#!/usr/bin/env python3
"""Create a bounded, lossless-for-included-fields public trajectory release.

Only completed registered runs are read. Personal machine provenance and
unrelated diagnostic channels are excluded; original files are untouched.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from quality_control import evaluate_runtime_diagnostics

FRAME_FIELDS = frozenset({
    "positions", "types", "ids", "radii", "core_radii", "physical_time_ms",
    "simulation_step", "delta_t_ms", "frame_interval_ms", "frame_write_complete",
    "total_force_magnitude", "total_force_vector", "electrostatic_force_vector",
    "hydrophobic_force_vector", "van_der_waals_force_vector",
    "steric_core_force_vector", "drift_scale", "drift_limited_particle_count",
    "previous_step_overlap_projection_corrected_particle_count",
    "previous_step_overlap_projection_iterations",
    "previous_step_overlap_projection_max_excess_nm",
    "previous_step_overlap_projection_pair_count",
    "previous_step_overlap_projection_residual_max_excess_nm",
    "previous_step_overlap_projection_residual_pair_count",
    "previous_step_overlap_projection_vesicle_center_cutoff_nm",
    "previous_step_overlap_projection_vesicle_neighbor_cell_reach",
})
OTHER_GROUPS = ("run_metadata", "run_completion", "minimizer_diagnostics")


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def copy_dataset(src: h5py.Dataset, target: h5py.Group, key: str) -> None:
    values = src[()]
    if src.shape == ():
        target.create_dataset(key, data=values)
    else:
        target.create_dataset(key, data=values, compression="gzip", compression_opts=6,
                              shuffle=True, chunks=True)


def export_run(source: Path, destination: Path, spec: dict) -> dict:
    expected_config = ROOT / spec["config"]
    if sha_file(expected_config) != spec["config_sha256"]:
        raise ValueError("Registered config checksum changed")
    if sha_file(source / "Dashboard.input.json") != spec["config_sha256"]:
        raise ValueError("Completed run used a different input config")
    if evaluate_runtime_diagnostics(source, 40000)["status"] != "passed":
        raise ValueError("Runtime diagnostics failed")
    candidates = list(source.glob("*.h5"))
    if len(candidates) != 1:
        raise ValueError("Expected exactly one HDF5 trajectory")
    original = candidates[0]
    destination.mkdir(parents=True, exist_ok=True)
    final_output = destination / "trajectory.h5"
    output = destination / "trajectory.h5.tmp"
    if final_output.exists() or output.exists():
        raise ValueError("Destination exists: refusing to overwrite reviewed data")
    with h5py.File(original, "r") as src, h5py.File(output, "x") as dst:
        frames = sorted((k for k in src if k.startswith("step_")), key=lambda k: int(k.split("_")[1]))
        if frames != [f"step_{step}" for step in range(200, 40001, 200)]:
            raise ValueError("Saved schedule is incomplete")
        if not np.asarray(src["run_completion/run_completed"][()]).item():
            raise ValueError("Run completion marker absent")
        for group_name in OTHER_GROUPS:
            original_group = src[group_name]
            new_group = dst.create_group(group_name)
            for key, dataset in original_group.items():
                if key.startswith(("sk_rock_", "srock2_")):
                    continue
                copy_dataset(dataset, new_group, key)
        for frame in frames:
            original_group = src[frame]
            missing = {"positions", "types", "ids", "radii", "physical_time_ms"} - set(original_group)
            if missing:
                raise ValueError("Missing required scientific frame fields")
            target_group = dst.create_group(frame)
            for key in sorted(FRAME_FIELDS & set(original_group)):
                copy_dataset(original_group[key], target_group, key)
        dst.attrs["publication_subset_schema"] = 1
        dst.attrs["publication_subset_scope"] = "registered force-response experiment only"
        dst.flush()
    # End-to-end exact equality of *every* copied dataset: no lossy resampling.
    with h5py.File(original, "r") as src, h5py.File(output, "r") as dst:
        def verify(name, object):
            if isinstance(object, h5py.Dataset):
                assert object.dtype == dst[name].dtype and object.shape == dst[name].shape
                if not np.array_equal(object[()], dst[name][()], equal_nan=np.issubdtype(object.dtype, np.floating)):
                    raise ValueError("Exported values differ: " + name)
        dst.visititems(lambda name, obj: verify(name, src[name]) if isinstance(obj, h5py.Dataset) else None)
    output.replace(final_output)
    shutil.copyfile(expected_config, destination / "Dashboard.input.json")
    diagnostics = json.loads((source / "numerical_diagnostics.json").read_text())
    (destination / "numerical_diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    return {"run_id": spec["run_id"], "seed": spec["seed"],
            "E": spec["electrostatic_enabled"], "H": spec["hydrophobic_enabled"],
            "status": "complete_runtime_qc_passed_values_verified",
            "original_h5_sha256": sha_file(original), "publication_h5_sha256": sha_file(final_output),
            "publication_h5_size_bytes": final_output.stat().st_size,
            "config_sha256": spec["config_sha256"], "saved_frames": 200}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="Private builds directory containing completed registered runs")
    parser.add_argument("--output", type=Path, default=ROOT / "data/raw")
    args = parser.parse_args()
    study = json.loads((ROOT / "configs/manifest.json").read_text())
    rows = []
    for spec in study["runs"]:
        candidates = sorted(args.source.glob(spec["run_id"] + "----*"))
        if len(candidates) != 1:
            raise ValueError("Expected one completed source folder for " + spec["run_id"])
        row = export_run(candidates[0], args.output / spec["run_id"], spec)
        rows.append(row)
        print("VERIFIED", spec["run_id"], row["publication_h5_size_bytes"], flush=True)
    status = {"schema_version": 1, "registered_runs": len(study["runs"]), "verified_released_runs": len(rows),
              "selection": "registered completed numerical-QC-passing runs only",
              "excluded": "original machine provenance, unneeded per-frame diagnostics, RNG checkpoint state",
              "runs": rows}
    (ROOT / "data/availability.json").write_text(json.dumps(status, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
