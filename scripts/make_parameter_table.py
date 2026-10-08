#!/usr/bin/env python3
"""Render a complete human-readable inventory from the registered inputs."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKIP_SECTIONS = {"article_production", "confirmatory_design"}  # echoed run labels, not physical inputs
SOURCES = {
    "environment.volume": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6573206/",
    "vesicles.mean_core_radius_nm": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5815350/",
    "vesicles.std_core_radius_nm": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5815350/",
    "vesicles.core_radius_min_nm": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2849067/",
    "vesicles.core_radius_max_nm": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2849067/",
    "vesicles.corona_thickness_nm": "https://pub.dzne.de/record/273920/files/DZNE-2024-01394.pdf?subformat=pdfa",
    "synapsinI.radius": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11628805/",
    "synapsinI.surface_potential_mV": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8424992/",
    "vesicles.surface_potential_mV": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2852195/",
    "parameters.solvent_relative_permittivity": "https://www.nist.gov/publications/static-dielectric-constant-liquid-water-between-274-and-418-k-near-saturated-vapor-0",
    "forces_parms.hydrophobic_decay_length_nm": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2115549/",
    "forces_parms.hydrophobic_energy_syn_syn_kbt": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2115549/",
    "forces_parms.hydrophobic_energy_syn_ves_kbt": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2115549/",
    "parameters.hamaker": "https://www.jstage.jst.go.jp/article/jos/67/11/67_ess18088/_html/-char/en",
}
NOT_DIRECT = {"synapsinI.radius", "synapsinI.surface_potential_mV",
              "forces_parms.hydrophobic_decay_length_nm",
              "forces_parms.hydrophobic_energy_syn_syn_kbt",
              "forces_parms.hydrophobic_energy_syn_ves_kbt"}


def flatten(value, prefix=""):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from flatten(item, prefix + ("." if prefix else "") + key)
    elif isinstance(value, list) and len(value) == 1 and isinstance(value[0], dict):
        yield from flatten(value[0], prefix)
    else:
        yield prefix, value


def main():
    manifest = json.loads((ROOT / "configs/manifest.json").read_text())
    arms = {}
    for row in manifest["runs"]:
        if row["seed"] != 9310002:
            continue
        key = f"E{int(row['electrostatic_enabled'])}H{int(row['hydrophobic_enabled'])}"
        data = json.loads((ROOT / row["config"]).read_text())
        arms[key] = {name: value for name, value in flatten(data)
                     if name.split(".", 1)[0] not in SKIP_SECTIONS}
    if set(arms) != {"E0H0", "E0H1", "E1H0", "E1H1"}:
        raise ValueError("Expected all four registered arms")
    output = ROOT / "docs/parameters.csv"
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["parameter", "unit_or_encoding", "E0H0", "E0H1", "E1H0", "E1H1",
                         "evidence_status", "source_url"])
        for name in sorted(arms["E0H0"]):
            values = [arms[arm].get(name) for arm in ("E0H0", "E0H1", "E1H0", "E1H1")]
            if any(value != values[0] for value in values) and name not in {
                "forces_parms.enable_electrostatic", "forces_parms.enable_hydrophobic",
                "total.address", "total.seed", "total.actual_seed"}:
                raise ValueError("Unregistered parameter changes between force arms: " + name)
            unit = next((suffix for suffix in ("_nm", "_ms", "_mV", "_kbt", "_molar")
                         if name.endswith(suffix)), "input/internal; see docs/UNITS.md")
            status = ("model_prior; source constrains mechanism/scale, not exact value" if name in NOT_DIRECT
                      else "literature-supported scale; see parameter sources" if name in SOURCES
                      else "registered model choice/derived/configuration; not independently measured")
            writer.writerow([name, unit, *(json.dumps(value, ensure_ascii=False) for value in values),
                             status, SOURCES.get(name, "")])
    print(f"PARAMETERS {output}")


if __name__ == "__main__":
    main()
