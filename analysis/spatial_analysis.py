#!/usr/bin/env python3
"""Scientific, backward-compatible spatial analysis of simulation HDF5 files.

This is stage 01 of the analysis pipeline.  It intentionally creates no
figures: its responsibility is to produce traceable, physically defined tables
that all later publication figures can share.

Main guarantees
---------------
* Existing trajectories are never modified and do not need to be rerun.
* A run keeps its full folder identity (including timestamp) plus a content
  fingerprint, while ``Condition_ID`` excludes the seed/run timestamp.
* One species-aware surface-gap graph is built per frame and reused for all
  cluster, composition, topology, and sensitivity metrics.
* Missing metadata are recorded explicitly; hidden scientific defaults are not
  silently presented as measured parameters.
* The simulation contract is strictly non-periodic: distances are direct
  Euclidean distances and no minimum-image or coordinate unwrapping is used.
* Output is partitioned by run and analysis configuration, so rerunning with a
  different threshold cannot overwrite a scientifically distinct result.

The compatibility columns ``time_step``, ``phi``, ``radius_of_gyration``,
``mean_velocity_S`` and ``mean_velocity_V`` are retained for older downstream
scripts.  Their precise definitions are written to Analysis_Metadata_v3.json.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import re
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    import h5py
except ImportError:  # pragma: no cover - handled by _require_dependencies
    h5py = None  # type: ignore[assignment]

import numpy as np
import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover - handled by _require_dependencies
    pa = None  # type: ignore[assignment]
    pq = None  # type: ignore[assignment]

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.stats import wasserstein_distance

sys.path.insert(0, str(Path(__file__).resolve().parent))
from quality_control import evaluate_runtime_diagnostics, realized_radius_metrics


CODE_VERSION = "4.2.0"
MANIFEST_SCHEMA_VERSION = "sv-run-manifest-4.0"
OUTPUT_SCHEMA_VERSION = "sv-spatial-tables-4.0"
BOUNDARY_CONDITION = "closed_reflecting_no_flux"
STEP_RE = re.compile(r"(\d+)$")


def saved_schedule_complete(
    saved_steps: Sequence[int], expected_steps: int, stride: int,
    last_time_ms: float, dt_ms: float,
) -> bool:
    """Check all planned saved states; the current solver omits step zero."""
    if stride <= 0 or expected_steps < 0 or not saved_steps:
        return False
    scheduled = list(range(stride, expected_steps + 1, stride))
    valid_steps = list(saved_steps) == scheduled or list(saved_steps) == [0, *scheduled]
    return bool(
        valid_steps and math.isfinite(last_time_ms) and math.isfinite(dt_ms)
        and math.isclose(last_time_ms, saved_steps[-1] * dt_ms,
                         rel_tol=1e-7, abs_tol=1e-8)
    )
TIMESTAMP_RE = re.compile(r"----(\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})$")


def detected_logical_cpu_count() -> int:
    """Return the CPUs available to this process, respecting Linux affinity."""

    affinity = getattr(os, "sched_getaffinity", None)
    if callable(affinity):
        try:
            return max(1, len(affinity(0)))
        except OSError:
            pass
    return max(1, int(os.cpu_count() or 1))


def recommended_worker_count(logical_cpus: Optional[int] = None) -> int:
    """Use all detected logical CPUs except two, with one worker minimum."""

    detected = detected_logical_cpu_count() if logical_cpus is None else logical_cpus
    return max(1, int(detected or 1) - 2)

DATASET_ALIASES: Mapping[str, Tuple[str, ...]] = {
    "positions": ("positions", "position", "coords", "coordinates"),
    "types": ("types", "type", "particle_types", "species"),
    "radii": ("radii", "radius", "particle_radii"),
    "ids": ("ids", "id", "particle_ids", "particle_id"),
    "physical_time_ms": ("physical_time_ms", "time_ms"),
    "force": (
        "total_force_magnitude",
        "force_magnitude",
        "total_force",
        "forces",
        "force",
    ),
}

FOLDER_PARAMETER_MAP = {
    "cutoff": "cut_off",
    "rep": "control_V_rep",
    "att": "control_V_att",
    "vr": "control_V_rep",
    "va": "control_V_att",
    "hrep": "control_H_rep",
    "hatt": "control_H_att",
    "hr": "control_H_rep",
    "ha": "control_H_att",
    "elect": "adjusted_electro",
    "eadj": "adjusted_electro",
    "e": "adjusted_electro",
    "vadj": "adjusted_vdw",
    "hadj": "adjusted_hydro",
    "seed": "seed",
    "shape": "requested_shape",
}

CONDITION_KEYS = (
    "ionic_strength_molar",
    "debye_length_nm",
    "hydrophobic_energy_syn_syn_kbt",
    "hydrophobic_energy_syn_ves_kbt",
    "hydrophobic_energy_ves_ves_kbt",
    "hydrophobic_decay_length_nm",
    "enable_electrostatic",
    "enable_hydrophobic",
    "enable_van_der_waals",
    "minimum_surface_separation_nm",
    "core_allowed_overlap_nm",
    "overlap_projection_residual_tolerance_nm",
    "overlap_projection_max_iterations",
    "dynamic_drift_limit_fraction_of_radius",
    "environment_shape",
    "boundary_condition",
    "environment_volume_nm3",
    "num_synapsin",
    "num_vesicle",
    "temperature",
    "eta",
    "drag_syn",
    "drag_ves",
    "D_syn",
    "D_ves",
    "synapsin_radius",
    "vesicle_mean_radius",
    "vesicle_std_radius",
    "synapsin_surface_potential_mV",
    "vesicle_surface_potential_mV",
    "vesicle_core_radius_min_nm",
    "vesicle_core_radius_max_nm",
    "vesicle_corona_thickness_nm",
    "box_length",
    "box_width",
    "box_height",
    "time_step_ms",
    "integrator",
)

MICRO_COLUMNS = [
    "schema_version",
    "analysis_config_id",
    "Run_ID",
    "Condition_ID",
    "seed",
    "seed_status",
    "frame_index",
    "simulation_step",
    "time_step",
    "time_ms",
    "time_source",
    "boundary_condition",
    "particle_id_source",
    "analysis_shell_nm",
    "cluster_id",
    "component_label",
    "cluster_size",
    "cluster_fraction",
    "n_vesicles",
    "n_synapsins",
    "gyration_ellipsoid_volume_nm3",
    "particle_volume_nm3",
    "ellipsoid_occupancy_proxy",
    "packing_fraction_phi",
    "phi",
    "phi_proxy_in_unit_interval",
    "radius_of_gyration_nm",
    "radius_of_gyration",
    "shape_anisotropy",
    "core_shell_index",
    "core_shell_strength",
    "core_species",
    "type_assortativity",
    "segregation_index",
    "excess_segregation_index",
    "observed_cross_edge_fraction",
    "random_mixing_cross_fraction",
    "cross_contact_enrichment",
    "homophily_S",
    "homophily_V",
    "mean_neighbors_SS",
    "mean_neighbors_VV",
    "mean_neighbors_cross_S",
    "mean_neighbors_cross_V",
    "mean_force_S",
    "mean_force_V",
    "deterministic_drift_proxy_S",
    "deterministic_drift_proxy_V",
    "force_over_drag_S",
    "force_over_drag_V",
    "mean_velocity_S",
    "mean_velocity_V",
    "member_ids",
]

FRAME_COLUMNS = [
    "schema_version",
    "analysis_config_id",
    "Run_ID",
    "Condition_ID",
    "seed",
    "seed_status",
    "frame_index",
    "simulation_step",
    "time_step",
    "time_ms",
    "time_source",
    "boundary_condition",
    "particle_id_source",
    "analysis_shell_nm",
    "n_particles",
    "n_synapsins",
    "n_vesicles",
    "n_components_all",
    "n_valid_clusters",
    "valid_clustered_fraction",
    "largest_component_size",
    "largest_component_fraction",
    "monomer_fraction",
    "weight_average_finite_cluster_size",
    "vesicle_recruited_fraction",
    "synapsin_recruited_fraction",
    "largest_component_vesicle_fraction",
    "mean_valid_cluster_size",
    "adjacency_edges_SS",
    "adjacency_edges_SV",
    "adjacency_edges_VV",
    "mean_adjacency_degree",
    "overlap_fraction_among_adjacency_edges",
    "overlap_pair_count",
    "max_outer_overlap_depth_nm",
    "projection_residual_tolerance_nm",
    "overlap_pairs_above_allowed_plus_tolerance",
    "max_pair_allowance_excess_nm",
    "max_excess_beyond_projection_tolerance_nm",
]

SENSITIVITY_COLUMNS = [
    "schema_version",
    "analysis_config_id",
    "Run_ID",
    "Condition_ID",
    "seed",
    "seed_status",
    "frame_index",
    "simulation_step",
    "time_step",
    "time_ms",
    "analysis_shell_nm",
    "n_components_all",
    "monomer_fraction",
    "largest_component_size",
    "largest_component_fraction",
    "weight_average_finite_cluster_size",
    "adjacency_edges_SS",
    "adjacency_edges_SV",
    "adjacency_edges_VV",
    "largest_component_vesicle_fraction",
    "vesicle_recruited_fraction",
    "synapsin_recruited_fraction",
]

FEATURE_DEFINITIONS: Mapping[str, str] = {
    "structural_adjacency": "outer-surface gap g_ij = ||x_i-x_j||-R_i_outer-R_j_outer; an edge exists when g_ij <= analysis_shell_nm",
    "primary_analysis_shell": "2 nm; shell sensitivity is reported at 0, 2, 4, and 6 nm; these are graph definitions, not force cutoffs",
    "largest_component_fraction": "size of the largest connected component divided by all particles in the saved frame",
    "valid_clustered_fraction": "fraction of all particles in components with size >= min_cluster_size",
    "radius_of_gyration_nm": "sqrt(trace of the unweighted particle-center gyration tensor) using direct non-periodic coordinates",
    "shape_anisotropy": "relative shape anisotropy kappa^2 in [0,1] from gyration eigenvalues",
    "ellipsoid_occupancy_proxy": "sum of particle sphere volumes divided by the uniform-gyration ellipsoid volume; a morphology proxy, not a direct Voronoi/local packing fraction",
    "phi": "legacy alias of ellipsoid_occupancy_proxy; later figures must label it as a proxy",
    "core_shell_index": "(mean radial position of synapsin - mean radial position of vesicle) / Rg; positive means vesicle-enriched core",
    "core_shell_strength": "1-Wasserstein distance between species radial distributions divided by Rg",
    "type_assortativity": "Newman assortativity of the undirected binary-species contact graph; zero is composition/degree-matched mixing",
    "segregation_index": "compatibility alias of type_assortativity",
    "excess_segregation_index": "1 - observed cross-edge fraction / random-label cross-edge expectation at fixed composition",
    "deterministic_drift_proxy": "force magnitude divided by Stokes drag using recorded viscosity, radius and species drag factor; excludes Brownian displacement",
    "late_window": "last configured fraction of saved frames; never automatically called equilibrium",
    "weight_average_finite_cluster_size": "sum(n_c^2)/sum(n_c) after excluding one largest component; monomers remain included",
    "vesicle_recruited_fraction": "fraction of all vesicles belonging to a connected component of size at least two at the primary shell",
    "pair_overlap_allowance": "core_allowed_overlap_nm plus one vesicle-corona thickness for each vesicle in the pair; evaluated on outer radii",
    "max_pair_allowance_excess_nm": "maximum over primary-graph pairs of max(0, outer-overlap depth minus the pair-specific overlap allowance)",
    "max_excess_beyond_projection_tolerance_nm": "max(0, pair-allowance excess minus the recorded projection residual tolerance), with a 1e-6 nm analysis slack matching the solver numerical tolerance",
}


class IntakeError(RuntimeError):
    """Raw data or metadata are insufficient for the requested analysis."""


def _require_dependencies() -> None:
    missing: List[str] = []
    if h5py is None:
        missing.append("h5py")
    if pa is None or pq is None:
        missing.append("pyarrow")
    if missing:
        joined = ", ".join(missing)
        raise RuntimeError(
            f"Missing required package(s): {joined}. Install them in the pipeline "
            "environment before running stage 01."
        )


def _normalized_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def _first_mapping(value: Any) -> Dict[str, Any]:
    if isinstance(value, list):
        value = value[0] if value else {}
    return dict(value) if isinstance(value, dict) else {}


def _section(data: Mapping[str, Any], name: str) -> Dict[str, Any]:
    normalized = {_normalized_key(str(key)): value for key, value in data.items()}
    return _first_mapping(normalized.get(_normalized_key(name), {}))


def _mapping_value(mapping: Mapping[str, Any], *names: str) -> Any:
    normalized = {_normalized_key(str(key)): value for key, value in mapping.items()}
    for name in names:
        key = _normalized_key(name)
        if key in normalized:
            return normalized[key]
    return None


def _finite_number(value: Any) -> Optional[float]:
    if isinstance(value, (list, tuple, np.ndarray)) and np.size(value) == 1:
        value = np.asarray(value).reshape(-1)[0]
    if isinstance(value, (bool, np.bool_)) or value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _string_value(value: Any) -> Optional[str]:
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, (list, tuple, np.ndarray)) and np.size(value) == 1:
        value = np.asarray(value).reshape(-1)[0]
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _coerce_bool(value: Any) -> Optional[bool]:
    if isinstance(value, (list, tuple, np.ndarray)) and np.size(value) == 1:
        value = np.asarray(value).reshape(-1)[0]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and int(value) in {0, 1}:
        return bool(value)
    text = _string_value(value)
    if text is None:
        return None
    normalized = _normalized_key(text)
    if normalized in {"true", "yes", "on", "1"}:
        return True
    if normalized in {"false", "no", "off", "0"}:
        return False
    return None


def _step_number(name: str) -> int:
    match = STEP_RE.search(name)
    if not match:
        raise IntakeError(f"HDF5 group has no numeric step suffix: {name!r}")
    return int(match.group(1))


def _atomic_json_dump(payload: Mapping[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".json", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def discover_h5(run_folder: Path, preferred: Optional[str] = None) -> Path:
    if preferred:
        candidate = Path(preferred).expanduser()
        if not candidate.is_absolute():
            candidate = run_folder / candidate
        if not candidate.is_file():
            raise IntakeError(f"Requested HDF5 file does not exist: {candidate}")
        return candidate.resolve()
    candidates = sorted(
        [*run_folder.glob("*.h5"), *run_folder.glob("*.hdf5")],
        key=lambda path: ("simulation_Rawdata" not in path.name, -path.stat().st_size, path.name),
    )
    if not candidates:
        raise IntakeError(f"No .h5/.hdf5 file found in {run_folder}")
    return candidates[0].resolve()


def resolve_dataset_names(group: Any) -> Dict[str, Optional[str]]:
    available = {_normalized_key(name): name for name in group.keys()}
    return {
        logical_name: next(
            (
                available[_normalized_key(alias)]
                for alias in aliases
                if _normalized_key(alias) in available
            ),
            None,
        )
        for logical_name, aliases in DATASET_ALIASES.items()
    }


def _valid_force_shape(shape: Tuple[int, ...], count: int) -> bool:
    return shape in {(count,), (count, 1), (count, 3)}


def _root_attribute_parameters(h5: Any) -> Dict[str, Any]:
    attrs = {_normalized_key(str(key)): value for key, value in h5.attrs.items()}

    def attr(*names: str) -> Any:
        for name in names:
            key = _normalized_key(name)
            if key in attrs:
                return attrs[key]
        return None

    return {
        "time_step_ms": _finite_number(attr("time_step_ms", "dt_ms")),
        "box_length": _finite_number(attr("box_length", "length", "lx")),
        "box_width": _finite_number(attr("box_width", "width", "ly")),
        "box_height": _finite_number(attr("box_height", "height", "lz")),
        "connectivity_SS_nm": _finite_number(attr("connectivity_SS_nm", "synapsin_synapsin")),
        "connectivity_SV_nm": _finite_number(attr("connectivity_SV_nm", "synapsin_vesicle")),
        "connectivity_VV_nm": _finite_number(attr("connectivity_VV_nm", "vesicle_vesicle")),
    }


def validate_hdf5(h5_path: Path) -> Dict[str, Any]:
    warnings: List[str] = []
    digest = hashlib.sha256()
    digest.update(h5_path.name.encode("utf-8"))
    digest.update(str(h5_path.stat().st_size).encode("ascii"))
    with h5py.File(h5_path, "r") as h5:  # type: ignore[union-attr]
        step_names = sorted(
            (name for name in h5.keys() if str(name).startswith("step_")),
            key=_step_number,
        )
        if not step_names:
            raise IntakeError("HDF5 contains no top-level step_* groups")
        step_values = [_step_number(name) for name in step_names]
        if len(set(step_values)) != len(step_values):
            raise IntakeError("Duplicate numeric simulation steps were found")

        reference_count: Optional[int] = None
        reference_names: Optional[Dict[str, Optional[str]]] = None
        particle_counts: List[int] = []
        type_values: set[int] = set()
        first_ids: Optional[np.ndarray] = None
        last_ids: Optional[np.ndarray] = None

        for frame_index, step_name in enumerate(step_names):
            group = h5[step_name]
            names = resolve_dataset_names(group)
            missing = [key for key in ("positions", "types", "radii") if names[key] is None]
            if missing:
                raise IntakeError(f"{step_name} is missing required datasets: {', '.join(missing)}")
            positions = group[names["positions"]]
            types = group[names["types"]]
            radii = group[names["radii"]]
            if len(positions.shape) != 2 or positions.shape[1] != 3:
                raise IntakeError(f"{step_name}/positions must have shape (N, 3), got {positions.shape}")
            count = int(positions.shape[0])
            if count == 0:
                raise IntakeError(f"{step_name} contains zero particles")
            if int(np.prod(types.shape)) != count or int(np.prod(radii.shape)) != count:
                raise IntakeError(f"{step_name} has inconsistent positions/types/radii lengths")
            if names["ids"] is not None and int(np.prod(group[names["ids"]].shape)) != count:
                raise IntakeError(f"{step_name} has an IDs array with the wrong length")
            if names["force"] is not None and not _valid_force_shape(group[names["force"]].shape, count):
                raise IntakeError(
                    f"{step_name}/{names['force']} must have shape (N,), (N,1), or (N,3); "
                    f"got {group[names['force']].shape}"
                )

            if reference_count is None:
                reference_count = count
            elif count != reference_count:
                warnings.append("particle_count_changes_between_frames")
            if reference_names is None:
                reference_names = names
            elif names != reference_names:
                warnings.append("dataset_names_change_between_frames")
            particle_counts.append(count)

            sample_n = min(64, count)
            pos_sample = np.asarray(positions[:sample_n], dtype=np.float64)
            rad_sample = np.asarray(radii[:sample_n], dtype=np.float64).reshape(-1)
            type_sample = np.asarray(types[:sample_n], dtype=np.int64).reshape(-1)
            if not np.all(np.isfinite(pos_sample)) or not np.all(np.isfinite(rad_sample)):
                raise IntakeError(f"{step_name} contains non-finite sampled positions/radii")
            if np.any(rad_sample <= 0):
                raise IntakeError(f"{step_name} contains non-positive sampled radii")
            type_values.update(int(value) for value in np.unique(type_sample))
            digest.update(step_name.encode("utf-8"))
            digest.update(pos_sample.tobytes())
            digest.update(rad_sample.tobytes())
            digest.update(type_sample.tobytes())

            if names["ids"] is not None and frame_index in {0, len(step_names) - 1}:
                ids = np.asarray(group[names["ids"]][:], dtype=np.int64).reshape(-1)
                if len(np.unique(ids)) != count:
                    raise IntakeError(f"{step_name} contains duplicate particle IDs")
                if frame_index == 0:
                    first_ids = ids
                else:
                    last_ids = ids

        identity_quality = (
            "recorded_particle_ids"
            if reference_names is not None and reference_names["ids"] is not None
            else "row_order_surrogate_ids"
        )
        ids_stable: Optional[bool] = None
        if identity_quality == "recorded_particle_ids":
            ids_stable = first_ids is not None and last_ids is not None and np.array_equal(first_ids, last_ids)
            if not ids_stable:
                warnings.append("particle_ids_not_stable_between_first_and_last_frame")
        else:
            warnings.append("particle_ids_missing_lineage_requires_verified_row_order")
        if reference_names is not None and reference_names["force"] is None:
            warnings.append("force_dataset_missing_drift_metrics_will_be_nan")

        h5_parameters = _root_attribute_parameters(h5)

    step_deltas = np.diff(step_values)
    stride = float(np.median(step_deltas)) if len(step_deltas) else None
    if len(step_deltas) and not np.allclose(step_deltas, stride):
        warnings.append("irregular_saved_frame_stride")
    return {
        "h5_filename": h5_path.name,
        "h5_size_bytes": int(h5_path.stat().st_size),
        "frame_count": len(step_names),
        "first_step": step_values[0],
        "last_step": step_values[-1],
        "saved_stride_steps": stride,
        "particle_count_min": min(particle_counts),
        "particle_count_max": max(particle_counts),
        "particle_types_sampled": sorted(type_values),
        "dataset_names": reference_names,
        "particle_identity_quality": identity_quality,
        "particle_ids_stable_first_to_last": ids_stable,
        "root_attribute_parameters": h5_parameters,
        "quick_content_fingerprint": digest.hexdigest(),
        "warnings": sorted(set(warnings)),
    }


def inspect_hdf5_light(h5_path: Path) -> Dict[str, Any]:
    """Read only the frame index plus first/last metadata needed by analysis.

    SimulationSV writes a fixed schema and the scientific arrays are consumed
    again during metric calculation.  Re-reading a sample from every frame at
    intake would duplicate I/O without changing a result, so the one-click
    pipeline performs this constant-cost inspection instead.
    """

    warnings: List[str] = []
    stat = h5_path.stat()
    with h5py.File(h5_path, "r") as h5:  # type: ignore[union-attr]
        step_names = sorted(
            (name for name in h5.keys() if str(name).startswith("step_")),
            key=_step_number,
        )
        if not step_names:
            raise IntakeError("HDF5 contains no top-level step_* groups")
        step_values = [_step_number(name) for name in step_names]
        first_group = h5[step_names[0]]
        last_group = h5[step_names[-1]]
        names = resolve_dataset_names(first_group)
        missing = [key for key in ("positions", "types", "radii") if names[key] is None]
        if missing:
            raise IntakeError(f"First saved frame is missing: {', '.join(missing)}")
        first_count = int(first_group[names["positions"]].shape[0])
        last_names = resolve_dataset_names(last_group)
        last_count = int(last_group[last_names["positions"]].shape[0])
        types = np.asarray(first_group[names["types"]][:], dtype=np.int64).reshape(-1)
        identity_quality = "recorded_particle_ids" if names["ids"] is not None else "row_order_surrogate_ids"
        ids_stable: Optional[bool] = None
        if names["ids"] is not None and last_names["ids"] is not None:
            first_ids = np.asarray(first_group[names["ids"]][:], dtype=np.int64).reshape(-1)
            last_ids = np.asarray(last_group[last_names["ids"]][:], dtype=np.int64).reshape(-1)
            ids_stable = bool(np.array_equal(first_ids, last_ids))
        step_deltas = np.diff(step_values)
        stride = float(np.median(step_deltas)) if len(step_deltas) else None
        digest = hashlib.sha256(
            json.dumps(
                {
                    "name": h5_path.name,
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                    "first_step": step_values[0],
                    "last_step": step_values[-1],
                    "first_count": first_count,
                    "last_count": last_count,
                },
                sort_keys=True,
            ).encode("utf-8")
        ).hexdigest()
        return {
            "h5_filename": h5_path.name,
            "h5_size_bytes": int(stat.st_size),
            "frame_count": len(step_names),
            "first_step": step_values[0],
            "last_step": step_values[-1],
            "saved_stride_steps": stride,
            "particle_count_min": min(first_count, last_count),
            "particle_count_max": max(first_count, last_count),
            "particle_types_sampled": sorted(int(value) for value in np.unique(types)),
            "dataset_names": names,
            "particle_identity_quality": identity_quality,
            "particle_ids_stable_first_to_last": ids_stable,
            "root_attribute_parameters": _root_attribute_parameters(h5),
            "quick_content_fingerprint": digest,
            "warnings": warnings,
        }


def parse_folder_parameters(folder_name: str) -> Dict[str, float]:
    scenario_name = TIMESTAMP_RE.sub("", folder_name)
    parsed: Dict[str, float] = {}
    keys = sorted(FOLDER_PARAMETER_MAP, key=len, reverse=True)
    pattern = re.compile(
        rf"(?:^|_)({'|'.join(re.escape(key) for key in keys)})-(-?(?:\d+(?:\.\d*)?|\.\d+))(?=_|$)",
        flags=re.IGNORECASE,
    )
    for match in pattern.finditer(scenario_name):
        parsed[FOLDER_PARAMETER_MAP[match.group(1).lower()]] = float(match.group(2))
    return parsed


def extract_dashboard_parameters(data: Mapping[str, Any]) -> Dict[str, Any]:
    forces = _section(data, "forces_parms")
    environment = _section(data, "environment")
    synapsin = _section(data, "synapsinI")
    vesicles = _section(data, "vesicles")
    parameters = _section(data, "parameters")
    times = _section(data, "times")
    total = _section(data, "total")
    primary_shell = _finite_number(_mapping_value(forces, "connectivity_shell_nm"))
    return {
        "adjusted_electro": _finite_number(_mapping_value(forces, "adjusted_electro")),
        "adjusted_hydro": _finite_number(_mapping_value(forces, "adjusted_hydro")),
        "adjusted_vdw": _finite_number(_mapping_value(forces, "adjusted_vdw")),
        "control_H_rep": _finite_number(_mapping_value(forces, "control_H_rep")),
        "control_H_att": _finite_number(_mapping_value(forces, "control_H_att")),
        "control_V_rep": _finite_number(_mapping_value(forces, "control_V_rep")),
        "control_V_att": _finite_number(_mapping_value(forces, "control_V_att")),
        "cut_off": _finite_number(_mapping_value(forces, "cut_off")),
        "analysis_surface_shell_nm": primary_shell,
        "hydrophobic_energy_syn_syn_kbt": _finite_number(
            _mapping_value(forces, "hydrophobic_energy_syn_syn_kbt")
        ),
        "hydrophobic_energy_syn_ves_kbt": _finite_number(
            _mapping_value(forces, "hydrophobic_energy_syn_ves_kbt")
        ),
        "hydrophobic_energy_ves_ves_kbt": _finite_number(
            _mapping_value(forces, "hydrophobic_energy_ves_ves_kbt")
        ),
        "hydrophobic_decay_length_nm": _finite_number(
            _mapping_value(forces, "hydrophobic_decay_length_nm")
        ),
        "enable_electrostatic": _coerce_bool(_mapping_value(forces, "enable_electrostatic")),
        "enable_hydrophobic": _coerce_bool(_mapping_value(forces, "enable_hydrophobic")),
        "enable_van_der_waals": _coerce_bool(_mapping_value(forces, "enable_van_der_waals")),
        "minimum_surface_separation_nm": _finite_number(
            _mapping_value(forces, "minimum_surface_separation_nm")
        ),
        "core_allowed_overlap_nm": _finite_number(
            _mapping_value(forces, "core_allowed_overlap_nm")
        ),
        "overlap_projection_residual_tolerance_nm": _finite_number(
            _mapping_value(forces, "overlap_projection_residual_tolerance_nm")
        ),
        "overlap_projection_max_iterations": _finite_number(
            _mapping_value(forces, "overlap_projection_max_iterations")
        ),
        "dynamic_drift_limit_fraction_of_radius": _finite_number(
            _mapping_value(forces, "dynamic_drift_limit_fraction_of_radius")
        ),
        "environment_shape": _finite_number(_mapping_value(environment, "shape")),
        "environment_volume_nm3": _finite_number(_mapping_value(environment, "volume")),
        "boundary_condition": _string_value(_mapping_value(environment, "boundary_condition")),
        "box_length": _finite_number(_mapping_value(environment, "length", "box_length", "lx")),
        "box_width": _finite_number(_mapping_value(environment, "width", "box_width", "ly")),
        "box_height": _finite_number(_mapping_value(environment, "height", "box_height", "lz")),
        "num_synapsin": _finite_number(_mapping_value(synapsin, "num")),
        "num_vesicle": _finite_number(_mapping_value(vesicles, "num")),
        "temperature": _finite_number(_mapping_value(parameters, "temperature")),
        "ionic_strength_molar": _finite_number(_mapping_value(parameters, "ionic_strength_molar")),
        "debye_length_nm": _finite_number(_mapping_value(parameters, "debye_length_nm")),
        "eta": _finite_number(_mapping_value(parameters, "eta")),
        "drag_syn": _finite_number(_mapping_value(synapsin, "drag_syn", "drag")),
        "drag_ves": _finite_number(_mapping_value(vesicles, "drag_ves", "drag")),
        "D_syn": _finite_number(_mapping_value(synapsin, "D_syn")),
        "D_ves": _finite_number(_mapping_value(vesicles, "D_ves")),
        "synapsin_radius": _finite_number(_mapping_value(synapsin, "radius")),
        "synapsin_surface_potential_mV": _finite_number(
            _mapping_value(synapsin, "surface_potential_mV")
        ),
        "vesicle_mean_radius": _finite_number(_mapping_value(vesicles, "mean_radius")),
        "vesicle_std_radius": _finite_number(_mapping_value(vesicles, "std_radius")),
        "vesicle_core_radius_min_nm": _finite_number(
            _mapping_value(vesicles, "core_radius_min_nm")
        ),
        "vesicle_core_radius_max_nm": _finite_number(
            _mapping_value(vesicles, "core_radius_max_nm")
        ),
        "vesicle_corona_thickness_nm": _finite_number(
            _mapping_value(vesicles, "corona_thickness_nm")
        ),
        "vesicle_surface_potential_mV": _finite_number(
            _mapping_value(vesicles, "surface_potential_mV")
        ),
        "synapsin_charge": _finite_number(_mapping_value(synapsin, "synapsinI_charge")),
        "vesicle_charge": _finite_number(_mapping_value(vesicles, "synaptophysin_charge")),
        "lambda_syn": _finite_number(_mapping_value(synapsin, "lambda_hydro")),
        "lambda_ves": _finite_number(_mapping_value(vesicles, "lambda_hydro")),
        "time_step_ms": _finite_number(_mapping_value(times, "time_step", "time_step_ms", "dt_ms")),
        "integrator": _string_value(_mapping_value(times, "integrator")),
        "run_label": _string_value(_mapping_value(total, "address")),
        # Kept as compatibility aliases.  The legacy 10/28/46 centre-distance
        # table is deliberately ignored; all three use the registered 2-nm
        # outer-surface structural shell.
        "connectivity_SS_nm": primary_shell,
        "connectivity_SV_nm": primary_shell,
        "connectivity_VV_nm": primary_shell,
    }


def _resolve_seed(data: Mapping[str, Any], folder_parameters: Mapping[str, float]) -> Dict[str, Any]:
    total = _section(data, "total")
    actual_seed = _finite_number(_mapping_value(total, "actual_seed"))
    configured_seed = _finite_number(_mapping_value(total, "seed"))
    raw_use_seed = _mapping_value(total, "use_seed")
    use_seed = _coerce_bool(raw_use_seed) is True
    if actual_seed is not None:
        return {"value": int(actual_seed), "status": "known", "source": "dashboard.actual_seed"}
    if use_seed and configured_seed is not None:
        return {"value": int(configured_seed), "status": "known", "source": "dashboard.configured_seed"}
    if "seed" in folder_parameters:
        return {"value": int(folder_parameters["seed"]), "status": "recovered", "source": "folder_name_fallback"}
    return {
        "value": None,
        "status": "unknown",
        "source": "not_recorded",
        "configured_seed_ignored": int(configured_seed) if configured_seed is not None else None,
        "reason": "actual runtime seed was not recorded; no artificial replicate identity was created",
    }


def build_run_manifest(
    run_folder: os.PathLike[str] | str,
    h5_filename: Optional[str] = None,
    dashboard_filename: str = "Dashboard.json",
    output_path: Optional[os.PathLike[str] | str] = None,
    parameter_overrides: Optional[Mapping[str, Any]] = None,
    require_dashboard: bool = False,
) -> Dict[str, Any]:
    _require_dependencies()
    run_path = Path(run_folder).expanduser().resolve()
    if not run_path.is_dir():
        raise IntakeError(f"Run folder does not exist: {run_path}")
    h5_path = discover_h5(run_path, h5_filename)
    dashboard_path = run_path / dashboard_filename
    dashboard: Dict[str, Any] = {}
    warnings: List[str] = []
    if dashboard_path.is_file():
        try:
            with dashboard_path.open("r", encoding="utf-8") as handle:
                loaded = json.load(handle)
            if not isinstance(loaded, dict):
                raise ValueError("Dashboard root is not a JSON object")
            dashboard = loaded
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            if require_dashboard:
                raise IntakeError(f"Cannot read valid Dashboard JSON: {dashboard_path}: {exc}") from exc
            warnings.append(f"dashboard_unreadable:{type(exc).__name__}")
    elif require_dashboard:
        raise IntakeError(f"Dashboard file does not exist: {dashboard_path}")
    else:
        warnings.append("dashboard_missing_using_hdf5_folder_and_explicit_metadata_only")

    h5_info = inspect_hdf5_light(h5_path)
    warnings.extend(h5_info.pop("warnings"))
    folder_name = run_path.name
    timestamp_match = TIMESTAMP_RE.search(folder_name)
    scenario_name = TIMESTAMP_RE.sub("", folder_name).strip()
    fingerprint = h5_info["quick_content_fingerprint"]
    run_id = f"{folder_name}--fp-{fingerprint[:12]}"

    dashboard_parameters = extract_dashboard_parameters(dashboard)
    h5_parameters = h5_info.get("root_attribute_parameters", {})
    effective_parameters = dict(dashboard_parameters)
    parameter_sources = {
        key: ("dashboard" if value is not None else "missing")
        for key, value in effective_parameters.items()
    }
    for key, value in h5_parameters.items():
        if value is not None and effective_parameters.get(key) is None:
            effective_parameters[key] = value
            parameter_sources[key] = "hdf5_root_attribute"

    folder_parameters = parse_folder_parameters(folder_name)
    conflicts: List[Dict[str, Any]] = []
    requested_shape = folder_parameters.get("requested_shape")
    for key, folder_value in folder_parameters.items():
        if key in {"seed", "requested_shape"}:
            continue
        previous = effective_parameters.get(key)
        if previous is not None and not math.isclose(float(previous), folder_value, rel_tol=1e-7, abs_tol=1e-9):
            conflicts.append(
                {
                    "parameter": key,
                    "previous_value": previous,
                    "folder_runtime_value": folder_value,
                    "selected_value": previous,
                    "reason": "Dashboard.input.json is authoritative; folder text is a label only",
                }
            )

    if requested_shape is not None:
        actual_shape = effective_parameters.get("environment_shape")
        if actual_shape is None or not math.isclose(float(actual_shape), requested_shape):
            warnings.append("shape_folder_label_not_verified_by_environment_metadata")

    requested_boundary = (parameter_overrides or {}).get("boundary_condition")
    if requested_boundary not in {None, "auto", BOUNDARY_CONDITION, False, 0}:
        raise IntakeError(
            "Stage 01 is intentionally non-periodic for this simulation project; "
            "a periodic boundary override is not allowed"
        )

    for key, value in (parameter_overrides or {}).items():
        if key == "boundary_condition":
            continue
        if value is None:
            continue
        previous = effective_parameters.get(key)
        if previous is not None and previous != value:
            conflicts.append(
                {
                    "parameter": key,
                    "previous_value": previous,
                    "explicit_override": value,
                    "selected_value": value,
                    "reason": "explicit stage-01 command-line override",
                }
            )
        effective_parameters[key] = value
        parameter_sources[key] = "explicit_override"

    recorded_boundary = effective_parameters.get("boundary_condition")
    if recorded_boundary not in {None, BOUNDARY_CONDITION}:
        raise IntakeError(
            f"Unsupported boundary condition {recorded_boundary!r}; current article analysis "
            f"is defined for {BOUNDARY_CONDITION!r}"
        )
    effective_parameters["boundary_condition"] = BOUNDARY_CONDITION
    parameter_sources["boundary_condition"] = (
        "dashboard" if recorded_boundary == BOUNDARY_CONDITION else "fixed_project_contract"
    )

    for key in ("connectivity_SS_nm", "connectivity_SV_nm", "connectivity_VV_nm"):
        value = effective_parameters.get(key)
        if value is None:
            raise IntakeError(
                f"{key} is missing. Publication analysis requires forces_parms.connectivity_shell_nm"
            )
        elif float(value) < 0:
            raise IntakeError(f"{key} must be non-negative, got {value}")

    seed = _resolve_seed(dashboard, folder_parameters)
    if seed["status"] == "unknown":
        warnings.append("actual_random_seed_unknown")
    expected_syn = effective_parameters.get("num_synapsin")
    expected_ves = effective_parameters.get("num_vesicle")
    if expected_syn is not None and expected_ves is not None:
        expected_total = int(round(float(expected_syn) + float(expected_ves)))
        if h5_info["particle_count_min"] != expected_total or h5_info["particle_count_max"] != expected_total:
            warnings.append(
                f"metadata_particle_total_{expected_total}_does_not_match_hdf5_"
                f"{h5_info['particle_count_min']}_{h5_info['particle_count_max']}"
            )

    condition_payload = {key: effective_parameters.get(key) for key in CONDITION_KEYS}
    canonical = json.dumps(condition_payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    condition_id = "cond-" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]
    dynamics_supported = (
        h5_info["particle_identity_quality"] == "recorded_particle_ids"
        and h5_info["particle_ids_stable_first_to_last"] is True
        and effective_parameters.get("time_step_ms") is not None
    )
    manifest: Dict[str, Any] = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "code_version": CODE_VERSION,
        "run_id": run_id,
        "legacy_run_label": scenario_name,
        "run_fingerprint": fingerprint,
        "condition_id": condition_id,
        "condition_parameters": condition_payload,
        "scenario_name": scenario_name,
        "folder_timestamp": timestamp_match.group(1) if timestamp_match else None,
        "seed": seed,
        "paths": {
            "run_folder": str(run_path),
            "h5": str(h5_path),
            "dashboard": str(dashboard_path.resolve()) if dashboard_path.is_file() else None,
        },
        "hdf5": h5_info,
        "dashboard_parameters": dashboard_parameters,
        "folder_recovered_parameters": parse_folder_parameters(folder_name),
        "effective_parameters": effective_parameters,
        "parameter_sources": parameter_sources,
        "metadata_conflicts": conflicts,
        "analysis_capabilities": {
            "spatial_connectivity": True,
            "particle_resolved_dynamics": dynamics_supported,
            "independent_replicate_identity": seed["status"] == "known",
            "periodic_geometry": False,
            "boundary_model": BOUNDARY_CONDITION,
        },
        "claim_scope": {
            "primary_label": "effective-interaction clustering/condensation analysis",
            "supported": [
                "species-aware geometric connectivity",
                "cluster-size and composition statistics",
                "time-resolved condensation proxies from saved trajectories",
            ],
            "requires_later_validation": [
                "equilibration or stationarity",
                "dense/dilute coexistence for thermodynamic LLPS",
                "run-level uncertainty across independent replicates",
            ],
            "not_identifiable_from_these_data_alone": [
                "independent mechanistic contribution of non-identifiable force terms",
                "replicate uncertainty when the runtime seed is unknown",
            ],
        },
        "warnings": sorted(set(warnings)),
    }
    destination = Path(output_path) if output_path else run_path / "Run_Manifest_v3.json"
    if not destination.is_absolute():
        destination = run_path / destination
    destination = destination.resolve()
    _atomic_json_dump(manifest, destination)
    manifest["manifest_path"] = str(destination)
    return manifest


def _nanmean(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    return float(np.mean(finite)) if finite.size else math.nan


def _read_force(group: Any, dataset_name: Optional[str], count: int) -> np.ndarray:
    if dataset_name is None:
        return np.full(count, np.nan, dtype=np.float64)
    force = np.asarray(group[dataset_name][:], dtype=np.float64)
    if force.shape == (count, 3):
        force = np.linalg.norm(force, axis=1)
    elif force.shape in {(count,), (count, 1)}:
        force = force.reshape(-1)
    else:
        raise ValueError(f"Force dataset has unsupported shape {force.shape}; expected (N,), (N,1), or (N,3)")
    return force


def _drift_proxy(
    force_magnitude: np.ndarray,
    radii: np.ndarray,
    types: np.ndarray,
    eta: Optional[float],
    drag_syn: Optional[float],
    drag_ves: Optional[float],
    synapsin_type: int,
    vesicle_type: int,
) -> np.ndarray:
    if eta is None or eta <= 0:
        return np.full(len(radii), np.nan, dtype=np.float64)
    drag = np.full(len(radii), np.nan, dtype=np.float64)
    syn = types == synapsin_type
    ves = types == vesicle_type
    if drag_syn is not None and drag_syn > 0:
        drag[syn] = 6.0 * np.pi * eta * radii[syn] * drag_syn
    if drag_ves is not None and drag_ves > 0:
        drag[ves] = 6.0 * np.pi * eta * radii[ves] * drag_ves
    result = np.full(len(radii), np.nan, dtype=np.float64)
    valid = np.isfinite(force_magnitude) & np.isfinite(drag) & (drag > 0)
    result[valid] = force_magnitude[valid] / drag[valid]
    return result


def _frame_context(task: Mapping[str, Any]) -> Dict[str, Any]:
    simulation_step = int(task["simulation_step"])
    dt = task.get("time_step_ms")
    recorded_time = _finite_number(task.get("physical_time_ms"))
    return {
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "analysis_config_id": task["analysis_config_id"],
        "Run_ID": task["run_id"],
        "Condition_ID": task["condition_id"],
        "seed": task.get("seed"),
        "seed_status": task["seed_status"],
        "frame_index": int(task["frame_index"]),
        "simulation_step": simulation_step,
        "time_step": simulation_step,
        "time_ms": (
            recorded_time
            if recorded_time is not None
            else float(simulation_step * dt) if dt is not None else math.nan
        ),
        "time_source": (
            "hdf5_physical_time_ms"
            if recorded_time is not None
            else "simulation_step_x_time_step_ms" if dt is not None else "simulation_step_only"
        ),
        "boundary_condition": BOUNDARY_CONDITION,
        "particle_id_source": task["particle_id_source"],
    }


def _type_assortativity(e_ss: int, e_sv: int, e_vv: int) -> float:
    total = e_ss + e_sv + e_vv
    if total == 0:
        return math.nan
    a_syn = (2.0 * e_ss + e_sv) / (2.0 * total)
    a_ves = (2.0 * e_vv + e_sv) / (2.0 * total)
    expected_trace = a_syn * a_syn + a_ves * a_ves
    denominator = 1.0 - expected_trace
    if denominator <= np.finfo(float).eps:
        return math.nan
    observed_trace = (e_ss + e_vv) / total
    return float((observed_trace - expected_trace) / denominator)


def _connected_components_from_pairs(
    count: int, pair_i: np.ndarray, pair_j: np.ndarray
) -> Tuple[int, np.ndarray, np.ndarray]:
    if len(pair_i):
        rows = np.concatenate([pair_i, pair_j])
        cols = np.concatenate([pair_j, pair_i])
        adjacency = csr_matrix(
            (np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(count, count)
        )
    else:
        adjacency = csr_matrix((count, count), dtype=np.uint8)
    n_components, labels = connected_components(adjacency, directed=False, return_labels=True)
    sizes = np.bincount(labels, minlength=n_components)
    return int(n_components), labels, sizes


def _weight_average_finite_cluster_size(component_sizes: np.ndarray) -> float:
    """Susceptibility-like finite-cluster size with one largest component removed."""

    sizes = np.asarray(component_sizes, dtype=np.float64)
    if sizes.size <= 1:
        return 0.0
    finite = np.delete(sizes, int(np.argmax(sizes)))
    denominator = float(np.sum(finite))
    return float(np.sum(finite * finite) / denominator) if denominator > 0 else 0.0


def _analyze_frame_arrays(
    task: Mapping[str, Any],
    positions: np.ndarray,
    types: np.ndarray,
    radii: np.ndarray,
    particle_ids: np.ndarray,
    force_magnitude: np.ndarray,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any], List[Dict[str, Any]]]:
    positions = np.asarray(positions, dtype=np.float64)
    types = np.asarray(types, dtype=np.int64).reshape(-1)
    radii = np.asarray(radii, dtype=np.float64).reshape(-1)
    particle_ids = np.asarray(particle_ids, dtype=np.int64).reshape(-1)
    force_magnitude = np.asarray(force_magnitude, dtype=np.float64).reshape(-1)
    count = len(positions)
    if positions.shape != (count, 3) or any(len(array) != count for array in (types, radii, particle_ids, force_magnitude)):
        raise ValueError(f"Inconsistent arrays in {task['step_name']}")
    if not np.all(np.isfinite(positions)) or not np.all(np.isfinite(radii)) or np.any(radii <= 0):
        raise ValueError(f"Invalid positions/radii in {task['step_name']}")
    if len(np.unique(particle_ids)) != count:
        raise ValueError(f"Duplicate particle IDs in {task['step_name']}")

    synapsin_type = int(task["synapsin_type"])
    vesicle_type = int(task["vesicle_type"])
    observed_types = set(int(value) for value in np.unique(types))
    allowed_types = {synapsin_type, vesicle_type}
    if not observed_types.issubset(allowed_types):
        raise ValueError(
            f"Unexpected particle type codes {sorted(observed_types - allowed_types)} in {task['step_name']}; "
            f"configured synapsin={synapsin_type}, vesicle={vesicle_type}"
        )

    force_over_drag = _drift_proxy(
        force_magnitude,
        radii,
        types,
        task.get("eta"),
        task.get("drag_syn"),
        task.get("drag_ves"),
        synapsin_type,
        vesicle_type,
    )
    if task.get("boundary_condition", BOUNDARY_CONDITION) != BOUNDARY_CONDITION:
        raise ValueError("Only direct non-periodic geometry is valid for this project")
    tree = cKDTree(positions)
    shells_nm = sorted(float(value) for value in task["shells_nm"])
    primary_shell_nm = float(task["primary_shell_nm"])
    search_radius = 2.0 * float(np.max(radii)) + max(shells_nm)
    candidate_pairs = tree.query_pairs(search_radius, output_type="ndarray")
    if candidate_pairs.size:
        pair_i = candidate_pairs[:, 0].astype(np.int64, copy=False)
        pair_j = candidate_pairs[:, 1].astype(np.int64, copy=False)
        center_distance = np.linalg.norm(positions[pair_i] - positions[pair_j], axis=1)
    else:
        pair_i = np.empty(0, dtype=np.int64)
        pair_j = np.empty(0, dtype=np.int64)
        center_distance = np.empty(0, dtype=np.float64)
    pair_type_i = types[pair_i]
    pair_type_j = types[pair_j]
    pair_ss = (pair_type_i == synapsin_type) & (pair_type_j == synapsin_type)
    pair_vv = (pair_type_i == vesicle_type) & (pair_type_j == vesicle_type)
    pair_sv = ((pair_type_i == synapsin_type) & (pair_type_j == vesicle_type)) | (
        (pair_type_i == vesicle_type) & (pair_type_j == synapsin_type)
    )
    surface_distance = center_distance - radii[pair_i] - radii[pair_j]
    primary_edge = surface_distance <= primary_shell_nm
    edge_i = pair_i[primary_edge]
    edge_j = pair_j[primary_edge]
    edge_surface_distance = surface_distance[primary_edge]
    edge_type_i = types[edge_i]
    edge_type_j = types[edge_j]
    edge_ss = (edge_type_i == synapsin_type) & (edge_type_j == synapsin_type)
    edge_vv = (edge_type_i == vesicle_type) & (edge_type_j == vesicle_type)
    edge_sv = ((edge_type_i == synapsin_type) & (edge_type_j == vesicle_type)) | (
        (edge_type_i == vesicle_type) & (edge_type_j == synapsin_type)
    )
    overlap_depth = np.maximum(0.0, -edge_surface_distance)
    core_allowance = float(task["core_allowed_overlap_nm"])
    corona_thickness = float(task["vesicle_corona_thickness_nm"])
    pair_allowance = (
        core_allowance
        + corona_thickness * (edge_type_i == vesicle_type)
        + corona_thickness * (edge_type_j == vesicle_type)
    )
    allowance_excess = np.maximum(0.0, overlap_depth - pair_allowance)
    projection_tolerance = float(task["overlap_projection_residual_tolerance_nm"])
    analysis_slack_nm = 1e-6
    excess_beyond_tolerance = np.maximum(0.0, allowance_excess - projection_tolerance)
    excess_beyond_tolerance[excess_beyond_tolerance <= analysis_slack_nm] = 0.0

    n_components, component_labels, component_sizes = _connected_components_from_pairs(
        count, edge_i, edge_j
    )
    primary_min_size = int(task["min_cluster_size"])
    valid_components = np.flatnonzero(component_sizes >= primary_min_size)
    base = {**_frame_context(task), "analysis_shell_nm": primary_shell_nm}

    sensitivity_rows: List[Dict[str, Any]] = []
    for shell_nm in shells_nm:
        shell_edge = surface_distance <= shell_nm
        shell_i = pair_i[shell_edge]
        shell_j = pair_j[shell_edge]
        shell_n, shell_labels, shell_sizes = _connected_components_from_pairs(count, shell_i, shell_j)
        shell_type_i = types[shell_i]
        shell_type_j = types[shell_j]
        shell_ss = (shell_type_i == synapsin_type) & (shell_type_j == synapsin_type)
        shell_vv = (shell_type_i == vesicle_type) & (shell_type_j == vesicle_type)
        shell_sv = ~shell_ss & ~shell_vv
        shell_largest = int(shell_sizes.max()) if len(shell_sizes) else 0
        shell_recruited = shell_sizes[shell_labels] >= 2
        shell_largest_mask = shell_labels == int(np.argmax(shell_sizes))
        shell_syn = types == synapsin_type
        shell_ves = types == vesicle_type
        sensitivity_rows.append(
            {
                **{key: base[key] for key in (
                    "schema_version", "analysis_config_id", "Run_ID", "Condition_ID", "seed",
                    "seed_status", "frame_index", "simulation_step", "time_step", "time_ms"
                )},
                "analysis_shell_nm": float(shell_nm),
                "n_components_all": int(shell_n),
                "monomer_fraction": float(np.sum(shell_sizes == 1) / count),
                "largest_component_size": shell_largest,
                "largest_component_fraction": float(shell_largest / count),
                "largest_component_vesicle_fraction": float(np.sum(shell_largest_mask & shell_ves) / np.sum(shell_ves)) if np.any(shell_ves) else math.nan,
                "vesicle_recruited_fraction": float(np.sum(shell_recruited & shell_ves) / np.sum(shell_ves)) if np.any(shell_ves) else math.nan,
                "synapsin_recruited_fraction": float(np.sum(shell_recruited & shell_syn) / np.sum(shell_syn)) if np.any(shell_syn) else math.nan,
                "weight_average_finite_cluster_size": _weight_average_finite_cluster_size(shell_sizes),
                "adjacency_edges_SS": int(np.sum(shell_ss)),
                "adjacency_edges_SV": int(np.sum(shell_sv)),
                "adjacency_edges_VV": int(np.sum(shell_vv)),
            }
        )

    valid_particle_count = int(component_sizes[valid_components].sum()) if len(valid_components) else 0
    largest_size = int(component_sizes.max()) if len(component_sizes) else 0
    recruited_mask = component_sizes[component_labels] >= 2
    syn_all = types == synapsin_type
    ves_all = types == vesicle_type
    n_syn_all = int(np.sum(syn_all))
    n_ves_all = int(np.sum(ves_all))
    largest_label = int(np.argmax(component_sizes)) if len(component_sizes) else -1
    largest_mask = component_labels == largest_label
    frame_summary: Dict[str, Any] = {
        **base,
        "n_particles": count,
        "n_synapsins": n_syn_all,
        "n_vesicles": n_ves_all,
        "n_components_all": int(n_components),
        "n_valid_clusters": int(len(valid_components)),
        "valid_clustered_fraction": float(valid_particle_count / count),
        "largest_component_size": largest_size,
        "largest_component_fraction": float(largest_size / count),
        "monomer_fraction": float(np.sum(component_sizes == 1) / count),
        "weight_average_finite_cluster_size": _weight_average_finite_cluster_size(component_sizes),
        "vesicle_recruited_fraction": float(np.sum(recruited_mask & ves_all) / n_ves_all) if n_ves_all else math.nan,
        "synapsin_recruited_fraction": float(np.sum(recruited_mask & syn_all) / n_syn_all) if n_syn_all else math.nan,
        "largest_component_vesicle_fraction": float(np.sum(largest_mask & ves_all) / n_ves_all) if n_ves_all else math.nan,
        "mean_valid_cluster_size": float(np.mean(component_sizes[valid_components])) if len(valid_components) else 0.0,
        "adjacency_edges_SS": int(np.sum(edge_ss)),
        "adjacency_edges_SV": int(np.sum(edge_sv)),
        "adjacency_edges_VV": int(np.sum(edge_vv)),
        "mean_adjacency_degree": float(2.0 * len(edge_i) / count),
        "overlap_fraction_among_adjacency_edges": (
            float(np.mean(edge_surface_distance < 0.0)) if len(edge_surface_distance) else 0.0
        ),
        "overlap_pair_count": int(np.sum(overlap_depth > 0.0)),
        "max_outer_overlap_depth_nm": float(np.max(overlap_depth)) if len(overlap_depth) else 0.0,
        "projection_residual_tolerance_nm": projection_tolerance,
        "overlap_pairs_above_allowed_plus_tolerance": int(np.sum(excess_beyond_tolerance > 0.0)),
        "max_pair_allowance_excess_nm": float(np.max(allowance_excess)) if len(allowance_excess) else 0.0,
        "max_excess_beyond_projection_tolerance_nm": (
            float(np.max(excess_beyond_tolerance)) if len(excess_beyond_tolerance) else 0.0
        ),
    }

    clusters: List[Dict[str, Any]] = []
    for component in valid_components:
        member_index = np.flatnonzero(component_labels == component)
        member_mask = np.zeros(count, dtype=bool)
        member_mask[member_index] = True
        internal_edge = member_mask[edge_i] & member_mask[edge_j]
        e_ss = int(np.sum(internal_edge & edge_ss))
        e_vv = int(np.sum(internal_edge & edge_vv))
        e_sv = int(np.sum(internal_edge & edge_sv))
        e_total = e_ss + e_vv + e_sv

        c_positions = positions[member_index].copy()
        c_types = types[member_index]
        c_radii = radii[member_index]
        c_force = force_magnitude[member_index]
        c_drift = force_over_drag[member_index]
        syn_mask = c_types == synapsin_type
        ves_mask = c_types == vesicle_type
        n_syn = int(np.sum(syn_mask))
        n_ves = int(np.sum(ves_mask))
        cluster_size = len(member_index)

        center = np.mean(c_positions, axis=0)
        centered = c_positions - center
        gyration = centered.T @ centered / cluster_size
        eigenvalues = np.sort(np.clip(np.linalg.eigvalsh(gyration), 0.0, None))[::-1]
        l1, l2, l3 = eigenvalues
        rg_sq = float(np.sum(eigenvalues))
        radius_of_gyration = math.sqrt(rg_sq)
        shape_anisotropy = (
            1.0 - 3.0 * float(l1 * l2 + l2 * l3 + l3 * l1) / (rg_sq * rg_sq)
            if rg_sq > 0
            else 0.0
        )
        ellipsoid_volume = float((4.0 * np.pi / 3.0) * math.sqrt(max(0.0, 125.0 * l1 * l2 * l3)))
        particle_volume = float(np.sum((4.0 * np.pi / 3.0) * c_radii**3))
        occupancy_proxy = particle_volume / ellipsoid_volume if ellipsoid_volume > 0 else math.nan
        phi_in_range = bool(0.0 <= occupancy_proxy <= 1.0) if math.isfinite(occupancy_proxy) else None

        radial_distance = np.linalg.norm(centered, axis=1)
        if n_syn and n_ves and radius_of_gyration > 0:
            mean_syn_radius = float(np.mean(radial_distance[syn_mask]))
            mean_ves_radius = float(np.mean(radial_distance[ves_mask]))
            core_shell_signed = (mean_syn_radius - mean_ves_radius) / radius_of_gyration
            core_shell_strength = float(
                wasserstein_distance(radial_distance[syn_mask], radial_distance[ves_mask])
                / radius_of_gyration
            )
            core_species = "vesicle" if core_shell_signed > 0 else "synapsin" if core_shell_signed < 0 else "none"
        else:
            core_shell_signed = math.nan
            core_shell_strength = math.nan
            core_species = "single_species"

        observed_cross = float(e_sv / e_total) if e_total else math.nan
        random_cross = (
            float(2.0 * n_syn * n_ves / (cluster_size * (cluster_size - 1)))
            if cluster_size > 1 and n_syn and n_ves
            else math.nan
        )
        cross_enrichment = (
            float(observed_cross / random_cross)
            if math.isfinite(observed_cross) and math.isfinite(random_cross) and random_cross > 0
            else math.nan
        )
        excess_segregation = 1.0 - cross_enrichment if math.isfinite(cross_enrichment) else math.nan
        assortativity = _type_assortativity(e_ss, e_sv, e_vv) if n_syn and n_ves else math.nan
        drift_syn = _nanmean(c_drift[syn_mask]) if n_syn else math.nan
        drift_ves = _nanmean(c_drift[ves_mask]) if n_ves else math.nan
        clusters.append(
            {
                **base,
                "cluster_id": f"{task['run_id']}::cfg_{task['analysis_config_id']}::frame_{task['frame_index']}::component_{int(component)}",
                "component_label": int(component),
                "cluster_size": cluster_size,
                "cluster_fraction": float(cluster_size / count),
                "n_vesicles": n_ves,
                "n_synapsins": n_syn,
                "gyration_ellipsoid_volume_nm3": ellipsoid_volume,
                "particle_volume_nm3": particle_volume,
                "ellipsoid_occupancy_proxy": occupancy_proxy,
                "packing_fraction_phi": occupancy_proxy,
                "phi": occupancy_proxy,
                "phi_proxy_in_unit_interval": phi_in_range,
                "radius_of_gyration_nm": radius_of_gyration,
                "radius_of_gyration": radius_of_gyration,
                "shape_anisotropy": float(shape_anisotropy),
                "core_shell_index": float(core_shell_signed),
                "core_shell_strength": core_shell_strength,
                "core_species": core_species,
                "type_assortativity": assortativity,
                "segregation_index": assortativity,
                "excess_segregation_index": excess_segregation,
                "observed_cross_edge_fraction": observed_cross,
                "random_mixing_cross_fraction": random_cross,
                "cross_contact_enrichment": cross_enrichment,
                "homophily_S": float((2 * e_ss) / (2 * e_ss + e_sv)) if (2 * e_ss + e_sv) else math.nan,
                "homophily_V": float((2 * e_vv) / (2 * e_vv + e_sv)) if (2 * e_vv + e_sv) else math.nan,
                "mean_neighbors_SS": float(2 * e_ss / n_syn) if n_syn else math.nan,
                "mean_neighbors_VV": float(2 * e_vv / n_ves) if n_ves else math.nan,
                "mean_neighbors_cross_S": float(e_sv / n_syn) if n_syn else math.nan,
                "mean_neighbors_cross_V": float(e_sv / n_ves) if n_ves else math.nan,
                "mean_force_S": _nanmean(c_force[syn_mask]) if n_syn else math.nan,
                "mean_force_V": _nanmean(c_force[ves_mask]) if n_ves else math.nan,
                "deterministic_drift_proxy_S": drift_syn,
                "deterministic_drift_proxy_V": drift_ves,
                "force_over_drag_S": drift_syn,
                "force_over_drag_V": drift_ves,
                "mean_velocity_S": drift_syn,
                "mean_velocity_V": drift_ves,
                "member_ids": particle_ids[member_index].tolist() if task["store_members"] else None,
            }
        )
    return clusters, frame_summary, sensitivity_rows


def _process_frame_from_open_h5(
    h5: Any, task: Mapping[str, Any]
) -> Tuple[List[Dict[str, Any]], Dict[str, Any], List[Dict[str, Any]]]:
    group = h5[task["step_name"]]
    names = resolve_dataset_names(group)
    missing = [key for key in ("positions", "types", "radii") if names[key] is None]
    if missing:
        raise ValueError(f"{task['step_name']} is missing datasets: {', '.join(missing)}")
    positions = np.asarray(group[names["positions"]][:], dtype=np.float64)
    types = np.asarray(group[names["types"]][:], dtype=np.int64).reshape(-1)
    radii = np.asarray(group[names["radii"]][:], dtype=np.float64).reshape(-1)
    count = len(positions)
    if names["ids"] is None:
        particle_ids = np.arange(1, count + 1, dtype=np.int64)
    else:
        particle_ids = np.asarray(group[names["ids"]][:], dtype=np.int64).reshape(-1)
    force_magnitude = _read_force(group, names["force"], count)
    frame_task = dict(task)
    if names.get("physical_time_ms") is not None:
        frame_task["physical_time_ms"] = _finite_number(
            np.asarray(group[names["physical_time_ms"]][()]).reshape(-1)[0]
        )
    return _analyze_frame_arrays(frame_task, positions, types, radii, particle_ids, force_magnitude)


def _process_task_batch(
    tasks: Sequence[Mapping[str, Any]],
) -> List[Tuple[List[Dict[str, Any]], Dict[str, Any], List[Dict[str, Any]]]]:
    if not tasks:
        return []
    with h5py.File(tasks[0]["h5_path"], "r") as h5:  # type: ignore[union-attr]
        return [_process_frame_from_open_h5(h5, task) for task in tasks]


def _chunked_tasks(tasks: Sequence[Mapping[str, Any]], workers: int) -> List[List[Mapping[str, Any]]]:
    if not tasks:
        return []
    batch_size = max(1, int(math.ceil(len(tasks) / max(1, workers * 2))))
    return [list(tasks[start : start + batch_size]) for start in range(0, len(tasks), batch_size)]


def _atomic_write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".parquet", dir=str(path.parent))
    os.close(fd)
    try:
        table = pa.Table.from_pandas(frame, preserve_index=False)  # type: ignore[union-attr]
        pq.write_table(table, temporary, compression="snappy")  # type: ignore[union-attr]
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _atomic_write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".csv", dir=str(path.parent))
    os.close(fd)
    try:
        frame.to_csv(temporary, index=False)
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _acquire_lock(lock_path: Path, timeout_seconds: float = 120.0) -> int:
    deadline = time.monotonic() + timeout_seconds
    while True:
        try:
            return os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            if time.monotonic() >= deadline:
                raise TimeoutError(f"Timed out waiting for output lock: {lock_path}")
            try:
                age = time.time() - lock_path.stat().st_mtime
                if age > max(600.0, timeout_seconds * 5.0):
                    lock_path.unlink()
                    continue
            except FileNotFoundError:
                continue
            time.sleep(0.1)


def _replace_analysis_in_csv(new_frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = path.with_suffix(path.suffix + ".lock")
    lock_fd = _acquire_lock(lock_path)
    try:
        if path.exists():
            existing = pd.read_csv(path)
            required = {"Run_ID", "analysis_config_id"}
            if not required.issubset(existing.columns):
                raise ValueError(f"Existing v3 summary lacks {sorted(required)} and will not be overwritten: {path}")
            keys = new_frame[["Run_ID", "analysis_config_id"]].astype(str).drop_duplicates()
            key_pairs = set(map(tuple, keys.to_numpy()))
            keep = [
                (str(run_id), str(config_id)) not in key_pairs
                for run_id, config_id in existing[["Run_ID", "analysis_config_id"]].to_numpy()
            ]
            combined = pd.concat([existing.loc[keep], new_frame], ignore_index=True, sort=False)
        else:
            combined = new_frame
        _atomic_write_csv(combined, path)
    finally:
        os.close(lock_fd)
        try:
            lock_path.unlink()
        except FileNotFoundError:
            pass


def _weighted_mean(frame: pd.DataFrame, value_column: str, weight_column: str) -> float:
    values = pd.to_numeric(frame[value_column], errors="coerce").to_numpy(dtype=float)
    weights = pd.to_numeric(frame[weight_column], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    return float(np.average(values[valid], weights=weights[valid])) if np.any(valid) else math.nan


class UnifiedSpatialEngine:
    """Extract the stage-01 spatial tables for one existing or future run."""

    def __init__(
        self,
        manifest: Mapping[str, Any],
        min_cluster_size: int = 2,
        workers: Optional[int] = None,
        store_members: bool = True,
        late_fraction: float = 0.30,
        primary_shell_nm: float = 2.0,
        shells_nm: Sequence[float] = (0.0, 2.0, 4.0, 6.0),
        synapsin_type: int = 1,
        vesicle_type: int = 2,
    ) -> None:
        if min_cluster_size < 2:
            raise ValueError("min_cluster_size must be at least 2")
        worker_limit = recommended_worker_count()
        requested_workers = worker_limit if workers is None else int(workers)
        if requested_workers < 1:
            raise ValueError("workers must be at least 1")
        if not 0 < late_fraction <= 1:
            raise ValueError("late_fraction must be in (0, 1]")
        if synapsin_type == vesicle_type:
            raise ValueError("synapsin_type and vesicle_type must be distinct")
        manifest_boundary = (
            manifest.get("effective_parameters", {}).get("boundary_condition")
            if isinstance(manifest.get("effective_parameters"), Mapping)
            else None
        )
        if manifest_boundary != BOUNDARY_CONDITION:
            raise ValueError(
                "The run manifest must declare "
                f"boundary_condition='{BOUNDARY_CONDITION}'; periodic geometry is outside "
                "this project's analysis contract"
            )
        self.manifest = dict(manifest)
        self.min_cluster_size = int(min_cluster_size)
        self.workers = min(requested_workers, worker_limit)
        self.store_members = bool(store_members)
        self.late_fraction = float(late_fraction)
        self.primary_shell_nm = float(primary_shell_nm)
        self.shells_nm = sorted({float(value) for value in shells_nm} | {self.primary_shell_nm})
        if self.primary_shell_nm < 0 or any(value < 0 for value in self.shells_nm):
            raise ValueError("analysis shells must be non-negative")
        registered_shell = _finite_number(
            manifest.get("effective_parameters", {}).get("analysis_surface_shell_nm")
        )
        if registered_shell is None or not math.isclose(
            registered_shell, self.primary_shell_nm, rel_tol=0.0, abs_tol=1e-12
        ):
            raise ValueError(
                "Primary analysis shell must equal forces_parms.connectivity_shell_nm "
                f"from Dashboard.input.json ({registered_shell!r} nm)"
            )
        self.synapsin_type = int(synapsin_type)
        self.vesicle_type = int(vesicle_type)
        config_payload = {
            "code_version": CODE_VERSION,
            "min_cluster_size": self.min_cluster_size,
            "store_members": self.store_members,
            "late_fraction": self.late_fraction,
            "synapsin_type": self.synapsin_type,
            "vesicle_type": self.vesicle_type,
            "boundary_condition": BOUNDARY_CONDITION,
            "primary_shell_nm": self.primary_shell_nm,
            "shells_nm": self.shells_nm,
        }
        canonical = json.dumps(config_payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
        self.analysis_config = config_payload
        self.analysis_config_id = "cfg-" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]
        self.df_micro = pd.DataFrame(columns=MICRO_COLUMNS)
        self.df_frames = pd.DataFrame(columns=FRAME_COLUMNS)
        self.df_sensitivity = pd.DataFrame(columns=SENSITIVITY_COLUMNS)
        self.df_macro = pd.DataFrame()

    def _tasks(self) -> List[Dict[str, Any]]:
        effective = self.manifest["effective_parameters"]
        h5_path = self.manifest["paths"]["h5"]
        with h5py.File(h5_path, "r") as h5:  # type: ignore[union-attr]
            step_names = sorted(
                (name for name in h5.keys() if str(name).startswith("step_")),
                key=_step_number,
            )
        common = {
            "h5_path": h5_path,
            "run_id": self.manifest["run_id"],
            "condition_id": self.manifest["condition_id"],
            "seed": self.manifest["seed"]["value"],
            "seed_status": self.manifest["seed"]["status"],
            "analysis_config_id": self.analysis_config_id,
            "time_step_ms": effective.get("time_step_ms"),
            "eta": effective.get("eta"),
            "drag_syn": effective.get("drag_syn"),
            "drag_ves": effective.get("drag_ves"),
            "core_allowed_overlap_nm": effective.get("core_allowed_overlap_nm"),
            "vesicle_corona_thickness_nm": effective.get("vesicle_corona_thickness_nm"),
            "overlap_projection_residual_tolerance_nm": effective.get(
                "overlap_projection_residual_tolerance_nm"
            ),
            "primary_shell_nm": self.primary_shell_nm,
            "shells_nm": self.shells_nm,
            "boundary_condition": BOUNDARY_CONDITION,
            "min_cluster_size": self.min_cluster_size,
            "store_members": self.store_members,
            "particle_id_source": self.manifest["hdf5"]["particle_identity_quality"],
            "synapsin_type": self.synapsin_type,
            "vesicle_type": self.vesicle_type,
        }
        if common["core_allowed_overlap_nm"] is None:
            raise ValueError("core_allowed_overlap_nm is required for pair-overlap diagnostics")
        if common["vesicle_corona_thickness_nm"] is None:
            raise ValueError("vesicle_corona_thickness_nm is required for pair-overlap diagnostics")
        if common["overlap_projection_residual_tolerance_nm"] is None:
            raise ValueError(
                "overlap_projection_residual_tolerance_nm is required for overlap diagnostics"
            )
        return [
            {
                **common,
                "step_name": step_name,
                "frame_index": frame_index,
                "simulation_step": _step_number(step_name),
            }
            for frame_index, step_name in enumerate(step_names)
        ]

    def execute(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        _require_dependencies()
        started = time.time()
        micro_rows: List[Dict[str, Any]] = []
        frame_rows: List[Dict[str, Any]] = []
        sensitivity_rows: List[Dict[str, Any]] = []
        tasks = self._tasks()
        print(
            f"[stage01] Run={self.manifest['run_id']} frames={len(tasks)} "
            f"config={self.analysis_config_id} workers={self.workers}"
        )
        batches = _chunked_tasks(tasks, self.workers)
        if self.workers == 1:
            batch_results: Iterable[Any] = map(_process_task_batch, batches)
            for batch in batch_results:
                for clusters, frame_summary, sensitivity in batch:
                    micro_rows.extend(clusters)
                    frame_rows.append(frame_summary)
                    sensitivity_rows.extend(sensitivity)
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.workers) as executor:
                for batch in executor.map(_process_task_batch, batches, chunksize=1):
                    for clusters, frame_summary, sensitivity in batch:
                        micro_rows.extend(clusters)
                        frame_rows.append(frame_summary)
                        sensitivity_rows.extend(sensitivity)

        self.df_micro = pd.DataFrame(micro_rows, columns=MICRO_COLUMNS)
        self.df_frames = pd.DataFrame(frame_rows, columns=FRAME_COLUMNS).sort_values("frame_index")
        self.df_sensitivity = pd.DataFrame(sensitivity_rows, columns=SENSITIVITY_COLUMNS).sort_values(
            ["frame_index", "analysis_shell_nm"]
        )
        self.df_macro = self._build_macro_summary()
        print(
            f"[stage01] Extracted {len(self.df_micro)} cluster rows, {len(self.df_frames)} frame rows, "
            f"and {len(self.df_sensitivity)} sensitivity rows in {time.time() - started:.1f} s"
        )
        return self.df_micro, self.df_frames, self.df_sensitivity, self.df_macro

    def _build_macro_summary(self) -> pd.DataFrame:
        frame_count = len(self.df_frames)
        late_count = max(1, int(math.ceil(frame_count * self.late_fraction)))
        late = self.df_frames.tail(late_count)
        effective = self.manifest["effective_parameters"]
        row: Dict[str, Any] = {
            "schema_version": OUTPUT_SCHEMA_VERSION,
            "code_version": CODE_VERSION,
            "analysis_config_id": self.analysis_config_id,
            "Run_ID": self.manifest["run_id"],
            "Legacy_Run_Label": self.manifest["legacy_run_label"],
            "Run_Fingerprint": self.manifest["run_fingerprint"],
            "Condition_ID": self.manifest["condition_id"],
            "seed": self.manifest["seed"]["value"],
            "seed_status": self.manifest["seed"]["status"],
            "scenario_name": self.manifest["scenario_name"],
            "frame_count": frame_count,
            "min_cluster_size": self.min_cluster_size,
            "primary_analysis_shell_nm": self.primary_shell_nm,
            "analysis_shells_nm": json.dumps(self.shells_nm),
            "late_window_fraction": self.late_fraction,
            "late_window_frame_count": late_count,
            "late_mean_largest_component_fraction": float(late["largest_component_fraction"].mean()),
            "late_mean_valid_clustered_fraction": float(late["valid_clustered_fraction"].mean()),
            "late_mean_monomer_fraction": float(late["monomer_fraction"].mean()),
            "late_mean_weight_average_finite_cluster_size": float(
                late["weight_average_finite_cluster_size"].mean()
            ),
            "late_mean_vesicle_recruited_fraction": float(late["vesicle_recruited_fraction"].mean()),
            "late_mean_synapsin_recruited_fraction": float(late["synapsin_recruited_fraction"].mean()),
            "late_mean_largest_component_vesicle_fraction": float(
                late["largest_component_vesicle_fraction"].mean()
            ),
            "late_mean_num_valid_clusters": float(late["n_valid_clusters"].mean()),
            "late_mean_overlap_fraction_among_adjacency_edges": float(
                late["overlap_fraction_among_adjacency_edges"].mean()
            ),
            "maximum_outer_overlap_depth_nm_over_saved_frames": float(
                self.df_frames["max_outer_overlap_depth_nm"].max()
            ),
            "maximum_pair_allowance_excess_nm_over_saved_frames": float(
                self.df_frames["max_pair_allowance_excess_nm"].max()
            ),
            "maximum_excess_beyond_projection_tolerance_nm_over_saved_frames": float(
                self.df_frames["max_excess_beyond_projection_tolerance_nm"].max()
            ),
            "frames_above_allowed_plus_tolerance_fraction": float(
                np.mean(self.df_frames["overlap_pairs_above_allowed_plus_tolerance"] > 0)
            ),
            "late_window_is_equilibrium": False,
            "equilibrium_status": "not_tested_late_window_only",
            "actual_seed_known": self.manifest["seed"]["status"] == "known",
            "metadata_warning_count": len(self.manifest.get("warnings", [])),
            "metadata_conflict_count": len(self.manifest.get("metadata_conflicts", [])),
            "claim_scope": self.manifest["claim_scope"]["primary_label"],
            "boundary_condition": effective["boundary_condition"],
            "phi_definition": "ellipsoid_occupancy_proxy_not_direct_packing_fraction",
            "segregation_definition": "newman_type_assortativity",
            "particle_member_ids_stored": self.store_members,
        }
        for key, value in effective.items():
            if isinstance(value, (str, int, float, bool)) or value is None:
                row[f"param_{key}"] = value

        dashboard_path = self.manifest["paths"].get("dashboard")
        expected_steps = None
        if dashboard_path:
            dashboard = json.loads(Path(dashboard_path).read_text(encoding="utf-8-sig"))
            total_steps = _mapping_value(_section(dashboard, "times"), "Ttime")
            if total_steps is not None:
                expected_steps = int(total_steps) - 1
        qc = evaluate_runtime_diagnostics(Path(self.manifest["paths"]["run_folder"]), expected_steps)
        row["runtime_qc_status"] = qc["status"]
        row["runtime_qc_reasons"] = json.dumps(qc["reasons"])
        for key, value in qc["observed"].items():
            row["runtime_" + key] = value
        row["saved_trajectory_complete"] = None
        if expected_steps is not None:
            saved_steps = self.df_frames["simulation_step"].astype(int).tolist()
            first_step = int(saved_steps[0])
            last_step = int(saved_steps[-1])
            last_time = float(self.df_frames["time_ms"].max())
            dt = effective.get("time_step_ms")
            stride = _mapping_value(_section(dashboard, "total"), "load")
            if stride is not None and int(stride) > 0 and dt is not None:
                # main.cpp saves fresh runs at load, 2*load, ... (not at step
                # zero).  Accept a future initial-frame export as well, but
                # require every scheduled saved step and no extra/missing one.
                row["saved_trajectory_complete"] = bool(
                    frame_count == len(saved_steps)
                    and saved_schedule_complete(
                        saved_steps, expected_steps, int(stride), last_time, float(dt)
                    )
                )
            row["saved_initial_state_present"] = first_step == 0
            row["saved_first_step"] = first_step
            row["saved_last_step"] = last_step
            row["saved_last_time_ms"] = last_time
        corona = effective.get("vesicle_corona_thickness_nm")
        volume = effective.get("environment_volume_nm3")
        if corona is not None and volume is not None:
            with h5py.File(self.manifest["paths"]["h5"], "r") as h5:
                first = h5[sorted((name for name in h5 if str(name).startswith("step_")), key=_step_number)[0]]
                names = resolve_dataset_names(first)
                row.update(realized_radius_metrics(first[names["radii"]][:], first[names["types"]][:], float(corona), float(volume)))

        if not self.df_micro.empty:
            late_frames = set(late["frame_index"].astype(int))
            late_micro = self.df_micro[self.df_micro["frame_index"].isin(late_frames)]
            row.update(
                {
                    "late_cluster_mean_ellipsoid_occupancy_proxy": float(late_micro["ellipsoid_occupancy_proxy"].mean()),
                    "late_particle_weighted_ellipsoid_occupancy_proxy": _weighted_mean(
                        late_micro, "ellipsoid_occupancy_proxy", "cluster_size"
                    ),
                    "late_particle_weighted_shape_anisotropy": _weighted_mean(
                        late_micro, "shape_anisotropy", "cluster_size"
                    ),
                    "late_particle_weighted_core_shell_signed": _weighted_mean(
                        late_micro, "core_shell_index", "cluster_size"
                    ),
                    "late_edge_structure_mean_type_assortativity": float(
                        late_micro["type_assortativity"].mean()
                    ),
                    "phi_proxy_out_of_range_cluster_count": int(
                        (late_micro["phi_proxy_in_unit_interval"] == False).sum()  # noqa: E712
                    ),
                }
            )
        return pd.DataFrame([row])

    def save(self, database_dir: os.PathLike[str] | str) -> Dict[str, Path]:
        database_path = Path(database_dir).expanduser().resolve()
        # Quick fingerprints are file-stat provenance, not unique run IDs.
        # Equal filenames/sizes/mtimes must not merge different trajectories.
        run_partition = hashlib.sha256((self.manifest["run_id"] + "|" + self.manifest["run_fingerprint"]).encode()).hexdigest()[:20]
        partition_path = database_path / "runs" / run_partition / self.analysis_config_id
        self.df_macro["run_partition"] = str(Path("runs") / run_partition / self.analysis_config_id)
        outputs = {
            "micro": partition_path / "Cluster_Timeseries_v4.parquet",
            "frames": partition_path / "Frame_Metrics_v4.parquet",
            "sensitivity": partition_path / "Shell_Sensitivity_v4.parquet",
            "run_summary": partition_path / "Run_Summary_v4.csv",
            "metadata": partition_path / "Analysis_Metadata_v4.json",
            "macro": database_path / "Run_Level_Metrics_v4.csv",
        }
        _atomic_write_parquet(self.df_micro, outputs["micro"])
        _atomic_write_parquet(self.df_frames, outputs["frames"])
        _atomic_write_parquet(self.df_sensitivity, outputs["sensitivity"])
        _atomic_write_csv(self.df_macro, outputs["run_summary"])
        metadata = {
            "schema_version": OUTPUT_SCHEMA_VERSION,
            "code_version": CODE_VERSION,
            "analysis_config_id": self.analysis_config_id,
            "analysis_config": self.analysis_config,
            "run_id": self.manifest["run_id"],
            "condition_id": self.manifest["condition_id"],
            "feature_definitions": dict(FEATURE_DEFINITIONS),
            "manifest_path": self.manifest.get("manifest_path"),
            "warnings": self.manifest.get("warnings", []),
            "metadata_conflicts": self.manifest.get("metadata_conflicts", []),
        }
        _atomic_json_dump(metadata, outputs["metadata"])
        _replace_analysis_in_csv(self.df_macro, outputs["macro"])
        for label, path in outputs.items():
            print(f"[stage01] {label}: {path}")
        return outputs


def _parse_shells(text: str) -> List[float]:
    try:
        values = sorted({float(part.strip()) for part in text.split(",") if part.strip()})
    except ValueError as exc:
        raise argparse.ArgumentTypeError("shells must be comma-separated numbers") from exc
    if not values or any((not math.isfinite(value)) or value < 0 for value in values):
        raise argparse.ArgumentTypeError("analysis shells must be finite and non-negative")
    return values


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_folder", help="Existing simulation output folder")
    parser.add_argument(
        "legacy_h5",
        nargs="?",
        help="Backward-compatible second positional HDF5 name; prefer --h5",
    )
    parser.add_argument("--h5", help="HDF5 filename or absolute path")
    parser.add_argument("--dashboard", default="Dashboard.input.json", help="authoritative input JSON filename")
    parser.add_argument("--require-dashboard", action="store_true")
    parser.add_argument("--manifest", help="Manifest path (default: RUN/Run_Manifest_v3.json)")
    parser.add_argument("--db-dir", help="Database directory (default: RUN/Analysis_v3)")
    parser.add_argument("--min-cluster-size", type=int, default=2)
    parser.add_argument("--primary-shell-nm", type=float, default=2.0)
    parser.add_argument("--shells-nm", type=_parse_shells, default=[0.0, 2.0, 4.0, 6.0])
    parser.add_argument("--workers", type=int, default=recommended_worker_count())
    parser.add_argument("--late-fraction", type=float, default=0.30)
    parser.add_argument("--no-member-ids", action="store_true", help="Reduce output size; disables overlap lineage later")
    parser.add_argument("--synapsin-type-code", type=int, default=1)
    parser.add_argument("--vesicle-type-code", type=int, default=2)
    parser.add_argument("--box-length", type=float)
    parser.add_argument("--box-width", type=float)
    parser.add_argument("--box-height", type=float)
    parser.add_argument("--time-step-ms", type=float)
    parser.add_argument("--connectivity-ss-nm", type=float)
    parser.add_argument("--connectivity-sv-nm", type=float)
    parser.add_argument("--connectivity-vv-nm", type=float)
    args = parser.parse_args(argv)

    if args.h5 and args.legacy_h5 and Path(args.h5).name != Path(args.legacy_h5).name:
        parser.error("Specify the HDF5 file either as the second positional argument or with --h5, not both")
    h5_filename = args.h5 or args.legacy_h5
    overrides = {
        "box_length": args.box_length,
        "box_width": args.box_width,
        "box_height": args.box_height,
        "time_step_ms": args.time_step_ms,
        "connectivity_SS_nm": args.connectivity_ss_nm,
        "connectivity_SV_nm": args.connectivity_sv_nm,
        "connectivity_VV_nm": args.connectivity_vv_nm,
    }
    run_folder = Path(args.run_folder).expanduser().resolve()
    manifest_path = Path(args.manifest).expanduser().resolve() if args.manifest else run_folder / "Run_Manifest_v3.json"
    database_dir = Path(args.db_dir).expanduser().resolve() if args.db_dir else run_folder / "Analysis_v3"
    try:
        manifest = build_run_manifest(
            run_folder,
            h5_filename=h5_filename,
            dashboard_filename=args.dashboard,
            output_path=manifest_path,
            parameter_overrides=overrides,
            require_dashboard=args.require_dashboard,
        )
        engine = UnifiedSpatialEngine(
            manifest,
            min_cluster_size=args.min_cluster_size,
            workers=args.workers,
            store_members=not args.no_member_ids,
            late_fraction=args.late_fraction,
            primary_shell_nm=args.primary_shell_nm,
            shells_nm=args.shells_nm,
            synapsin_type=args.synapsin_type_code,
            vesicle_type=args.vesicle_type_code,
        )
        engine.execute()
        engine.save(database_dir)
    except (IntakeError, RuntimeError, ValueError, OSError) as exc:
        print(f"[stage01:error] {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
