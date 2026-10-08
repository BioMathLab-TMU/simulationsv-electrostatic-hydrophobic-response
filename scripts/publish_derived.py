#!/usr/bin/env python3
"""Copy only verified, article-scoped results into the public data package."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FILES = (
    "analysis_status.json", "run_inventory.csv", "run_endpoints.csv",
    "frame_metrics.csv", "shell_sensitivity.csv", "condition_summary.csv",
    "paired_contrasts.csv",
)
FIGURES = (
    "assembly_growth", "registered_endpoints", "planned_paired_contrasts",
    "analysis_shell_sensitivity",
)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(4 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    status = json.loads((RESULTS / "analysis_status.json").read_text())
    available = json.loads((ROOT / "data/availability.json").read_text())
    if not (status["confirmatory_inference_ready"] and status["included_runs"] == 12
            and available["verified_released_runs"] == 12):
        raise ValueError("Do not publish incomplete/failing results")
    import csv
    with (RESULTS / "run_inventory.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 12 or any(row["status"] != "included" for row in rows):
        raise ValueError("Inventory and QC are not complete")
    processed = ROOT / "data/processed"
    processed.mkdir(parents=True, exist_ok=True)
    fig_target = ROOT / "figures"
    fig_target.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        shutil.copyfile(RESULTS / name, processed / name)
    for name in FIGURES:
        for extension in ("pdf", "png", "tif"):
            shutil.copyfile(RESULTS / "figures" / f"{name}.{extension}",
                            fig_target / f"{name}.{extension}")
    checks = {}
    for folder in (ROOT / "data/raw", processed, fig_target):
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                checks[path.relative_to(ROOT).as_posix()] = sha(path)
    (ROOT / "provenance/release_checksums.json").write_text(
        json.dumps({"schema_version": 1, "sha256": checks}, indent=2) + "\n")
    print(json.dumps({"raw_runs": 12, "processed_files": len(FILES),
                      "figure_files": len(FIGURES) * 3, "checksummed_files": len(checks)}))


if __name__ == "__main__":
    main()
