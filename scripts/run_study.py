#!/usr/bin/env python3
"""Run only the registered study sequentially; never alter its parameters."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from quality_control import evaluate_runtime_diagnostics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", default="build/SimulationSV")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--run-id", help="Run one registered configuration instead of the full sequence")
    args = parser.parse_args()
    if args.threads is not None and args.threads < 1:
        parser.error("threads must be positive")
    manifest = json.loads((ROOT / "configs/manifest.json").read_text())
    rows = [r for r in manifest["runs"] if args.run_id is None or args.run_id == r["run_id"]]
    if not rows:
        parser.error("run ID is not registered")
    # Verify the whole selected sequence before starting any computation.
    for row in rows:
        path = ROOT / row["config"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["config_sha256"]:
            raise ValueError("Configuration checksum mismatch: " + row["run_id"])
    if args.dry_run:
        print(json.dumps({"registered_runs": len(rows), "parameters_unchanged": True}))
        return 0
    executable = (ROOT / args.executable).resolve()
    if not executable.is_file():
        parser.error("executable does not exist; build it first")
    env = dict(os.environ, OMP_DYNAMIC="FALSE")
    if args.threads is not None:
        env["OMP_NUM_THREADS"] = str(args.threads)
    env["SIMULATIONSV_EXECUTABLE_SHA256"] = hashlib.sha256(executable.read_bytes()).hexdigest()
    for row in rows:
        existing = list((ROOT / "builds").glob(row["run_id"] + "----*"))
        if existing:
            if len(existing) != 1 or evaluate_runtime_diagnostics(existing[0], 40000)["status"] != "passed":
                raise RuntimeError("Existing run requires inspection; do not overwrite: " + row["run_id"])
            continue
        print("START " + row["run_id"], flush=True)
        subprocess.run([str(executable), str(ROOT / row["config"])], cwd=ROOT, env=env, check=True)
        folders = list((ROOT / "builds").glob(row["run_id"] + "----*"))
        if len(folders) != 1 or evaluate_runtime_diagnostics(folders[0], 40000)["status"] != "passed":
            raise RuntimeError("Recorded numerical QC failed: " + row["run_id"])
        print("COMPLETE " + row["run_id"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
