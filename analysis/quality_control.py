"""Read-only checks using the already approved article production guardrails.

These checks establish recorded numerical admissibility, not timestep
convergence, equilibration, or thermodynamic LLPS.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping


def check_diagnostics(d: Mapping[str, Any], expected_steps: int | None = None) -> dict[str, Any]:
    integer_fields = {
        "dynamic_steps_recorded", "nonfinite_euler_proposal_count",
        "nonfinite_position_update_count", "boundary_rejected_update_count",
        "projection_nonconverged_steps",
    }
    limits = {
        "maximum_projection_residual_excess_nm": 0.01 + 1e-12,
        "maximum_used_drift_fraction_of_radius": 0.2500001,
        "maximum_total_proposal_fraction_of_radius": 1.0,
        "maximum_brownian_displacement_sigma_units": 8.0,
    }
    required = integer_fields | set(limits)
    values: dict[str, float] = {}
    malformed = []
    for field in sorted(required):
        value = d.get(field)
        try:
            number = float(value)
            valid = not isinstance(value, bool) and math.isfinite(number) and number >= 0
            valid = valid and (field not in integer_fields or number.is_integer())
        except (TypeError, ValueError, OverflowError):
            valid = False
        if not valid:
            malformed.append(field)
        else:
            values[field] = number
    if malformed:
        return {"status": "failed", "reasons": ["missing_or_invalid:" + f for f in malformed], "checks": {}, "observed": values}
    steps = int(values["dynamic_steps_recorded"])
    checks = {field: values[field] <= limit for field, limit in limits.items()}
    checks.update({field: values[field] == 0 for field in integer_fields - {"dynamic_steps_recorded", "projection_nonconverged_steps"}})
    checks["dynamic_steps_positive"] = steps > 0
    if expected_steps is not None:
        checks["expected_dynamic_steps"] = steps == expected_steps
    checks["projection_nonconvergence_limit"] = values["projection_nonconverged_steps"] <= max(10, math.ceil(0.001 * steps))
    return {"status": "passed" if all(checks.values()) else "failed", "reasons": [key for key, passed in checks.items() if not passed], "checks": checks, "observed": values}


def evaluate_runtime_diagnostics(run_folder: Path, expected_steps: int | None = None) -> dict[str, Any]:
    path = Path(run_folder) / "numerical_diagnostics.json"
    if not path.is_file():
        return {"status": "unavailable", "reasons": ["numerical_diagnostics.json_missing"], "observed": {}, "checks": {}}
    try:
        data = json.loads(path.read_text(encoding="utf-8-sig"))
        if not isinstance(data, dict):
            raise ValueError("expected a JSON object")
        result = check_diagnostics(data, expected_steps)
    except (OSError, ValueError) as error:
        result = {"status": "failed", "reasons": ["unreadable_diagnostics:" + str(error)], "checks": {}, "observed": {}}
    result["source"] = str(path)
    return result


def realized_radius_metrics(outer_radii, types, corona_nm: float, volume_nm3: float) -> dict[str, Any]:
    import numpy as np

    radii = np.asarray(outer_radii, dtype=float)[np.asarray(types) == 2] - float(corona_nm)
    if not math.isfinite(volume_nm3) or volume_nm3 <= 0 or not math.isfinite(corona_nm) or corona_nm < 0:
        raise ValueError("Corona and compartment volume must be valid")
    if not np.isfinite(radii).all() or (radii <= 0).any():
        raise ValueError("Saved outer vesicle radii minus corona must be positive finite core radii")
    result: dict[str, Any] = {"realized_vesicle_count": int(len(radii)), "radius_reference": "first_saved_frame", "packing_definition": "nominal_sum_of_core_sphere_volumes_over_compartment_volume_not_union_volume"}
    for label, value in {
        "mean": np.mean(radii) if len(radii) else np.nan,
        "sd": np.std(radii, ddof=0) if len(radii) else np.nan,
        "min": np.min(radii) if len(radii) else np.nan,
        "max": np.max(radii) if len(radii) else np.nan,
        "third_moment": np.mean(radii**3) if len(radii) else np.nan,
    }.items():
        result["realized_vesicle_core_radius_" + label + ("_nm3" if label == "third_moment" else "_nm")] = float(value)
    result["realized_nominal_core_packing_fraction"] = float(4 * np.pi / 3 * np.sum(radii**3) / volume_nm3)
    return result
