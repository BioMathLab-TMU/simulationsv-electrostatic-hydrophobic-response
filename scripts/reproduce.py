#!/usr/bin/env python3
"""Reproduce registered spatial endpoints, paired contrasts, and figures."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from force_statistics import ENDPOINTS, paired_statistics
from plotting import COLORS, save_figure
from force_figures import make_figures
from spatial_analysis import UnifiedSpatialEngine, build_run_manifest
from quality_control import evaluate_runtime_diagnostics

STYLES = {(1, 1): (COLORS["both_inside"], "-"), (1, 0): (COLORS["ves_alone"], "--"),
          (0, 1): (COLORS["syn_alone"], "-."), (0, 0): (COLORS["neutral"], ":")}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "data/raw")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args(argv)
    args.input = args.input.expanduser().resolve()
    args.output = args.output.expanduser().resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    study = json.loads((ROOT / "configs/manifest.json").read_text())
    run_rows, frame_tables, shell_tables, inventory = [], [], [], []
    for spec in study["runs"]:
        folders = sorted(p for p in args.input.glob(spec["run_id"] + "*")
                         if p.is_dir() and (p.name == spec["run_id"] or p.name.startswith(spec["run_id"] + "----")))
        row = {"run_id": spec["run_id"], "seed": spec["seed"], "E": int(spec["electrostatic_enabled"]), "H": int(spec["hydrophobic_enabled"])}
        if not folders:
            inventory.append({**row, "status": "missing"})
            continue
        if len(folders) != 1:
            raise ValueError("Multiple copies of the same registered run: " + spec["run_id"])
        folder = folders[0]
        config = folder / "Dashboard.input.json"
        if not config.is_file() or hashlib.sha256(config.read_bytes()).hexdigest() != spec["config_sha256"]:
            inventory.append({**row, "status": "failed_config_integrity"})
            continue
        qc = evaluate_runtime_diagnostics(folder, 40000)
        if qc["status"] != "passed":
            inventory.append({**row, "status": "failed_runtime_qc", "reason": ";".join(qc["reasons"])})
            continue
        try:
            manifest = build_run_manifest(folder, dashboard_filename="Dashboard.input.json", require_dashboard=True,
                                          output_path=args.output / "cache" / spec["run_id"] / "manifest.json")
            engine = UnifiedSpatialEngine(manifest, workers=args.workers, primary_shell_nm=2,
                                          shells_nm=[0, 2, 4, 6], late_fraction=0.30, store_members=True)
            _, frames, shells, macro = engine.execute()
            # Keep the full extraction in the local reproducible cache.
            engine.save(args.output / "cache" / spec["run_id"])
            metrics = macro.iloc[0].to_dict()
            valid = bool(metrics.get("saved_trajectory_complete") and metrics.get("actual_seed_known"))
            valid = valid and int(manifest["seed"]["value"]) == int(spec["seed"])
            effective = manifest["effective_parameters"]
            valid = valid and effective.get("enable_electrostatic") == spec["electrostatic_enabled"]
            valid = valid and effective.get("enable_hydrophobic") == spec["hydrophobic_enabled"]
            valid = valid and effective.get("enable_van_der_waals") is True
            valid = valid and float(metrics["maximum_excess_beyond_projection_tolerance_nm_over_saved_frames"]) <= 0.01
            if not valid:
                inventory.append({**row, "status": "failed_saved_state_or_metadata_qc"})
                continue
            # Export publication tables without machine paths or private cache metadata.
            run_rows.append({**row, "electrostatic_enabled": row["E"], "hydrophobic_enabled": row["H"],
                             "van_der_waals_enabled": True, "frame_count": int(metrics["frame_count"]),
                             "late_window_frame_count": int(metrics["late_window_frame_count"]),
                             **{key: float(metrics[key]) for key in ENDPOINTS},
                             "mean_drift_limited_fraction": json.loads((folder / "numerical_diagnostics.json").read_text())["mean_drift_limited_fraction"]})
            for table, destination in ((frames, frame_tables), (shells, shell_tables)):
                selected = [c for c in table if c not in {"Run_ID", "Condition_ID", "analysis_config_id", "seed_status", "schema_version"}]
                destination.append(table[selected].assign(run_id=spec["run_id"], E=row["E"], H=row["H"]))
            inventory.append({**row, "status": "included"})
        except (ValueError, RuntimeError, OSError) as error:
            inventory.append({**row, "status": "failed_intake", "reason": type(error).__name__})
    pd.DataFrame(inventory).to_csv(args.output / "run_inventory.csv", index=False)
    ready = len(run_rows) == len(study["runs"]) == 12
    status = {"planned_runs": 12, "included_runs": len(run_rows), "confirmatory_inference_ready": ready,
              "independent_unit": "seed/run", "primary_shell_nm": 2, "late_fraction": 0.30,
              "scope": "finite-time model organization; not thermodynamic LLPS"}
    (args.output / "analysis_status.json").write_text(json.dumps(status, indent=2) + "\n")
    if run_rows:
        run_level = pd.DataFrame(run_rows)
        run_level.to_csv(args.output / "run_endpoints.csv", index=False)
        frames = pd.concat(frame_tables, ignore_index=True)
        shells = pd.concat(shell_tables, ignore_index=True)
        frames.to_csv(args.output / "frame_metrics.csv", index=False)
        shells.to_csv(args.output / "shell_sensitivity.csv", index=False)
        fig, ax = plt.subplots(figsize=(6, 4))
        for (e, h), table in frames.groupby(["E", "H"]):
            color, style = STYLES[(e, h)]
            for _, trajectory in table.groupby("seed"):
                ax.plot(trajectory["time_ms"], trajectory["largest_component_fraction"], color=color, alpha=0.25, lw=0.7)
            summary = table.groupby("time_ms")["largest_component_fraction"].mean()
            ax.plot(summary.index, summary.values, color=color, ls=style, label=f"E{e}H{h} (n={table['seed'].nunique()})")
        ax.set(xlabel="Model time (ms)", ylabel="Largest-component fraction (all particles)")
        ax.legend()
        save_figure(fig, args.output / "figures", "assembly_growth")
        if ready:
            summaries, contrasts = paired_statistics(run_level)
            summaries.to_csv(args.output / "condition_summary.csv", index=False)
            contrasts.to_csv(args.output / "paired_contrasts.csv", index=False)
            make_figures(run_level, contrasts, shells, args.output / "figures")
    print(json.dumps(status))
    return 0 if ready or not args.require_complete else 3


if __name__ == "__main__":
    raise SystemExit(main())
