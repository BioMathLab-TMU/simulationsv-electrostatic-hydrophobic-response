"""Mechanically generate S1 Text tables and a checksum-indexed S1 Data ZIP.

Reads article-specific results only. Manuscript output is written to the
canonical writing folder; no raw trajectories are copied into that workspace.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import zipfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def table(headers, rows, columns):
    return "\n".join([r"\begin{center}\small", r"\begin{tabular}{"+columns+"}",
        r"\toprule", " & ".join(headers)+r" \\ \midrule",
        *[" & ".join(row)+r" \\" for row in rows],
        r"\bottomrule\end{tabular}\end{center}"])


def atomic_text(path, content):
    fd, temporary = tempfile.mkstemp(prefix=".support_", suffix=path.suffix, dir=path.parent)
    os.close(fd)
    Path(temporary).write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manuscript", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    a.manuscript.mkdir(parents=True, exist_ok=True)
    endpoints = pd.read_csv(ROOT/"data/processed/run_endpoints.csv").sort_values(["E", "H", "seed"])
    interface = pd.read_csv(ROOT/"results/interface_diagnostics/interface_run_metrics.csv").sort_values(["E", "H", "seed"])
    mobility = pd.read_csv(ROOT/"results/mobility_diagnostic/run_level_mobility.csv")
    shells = pd.read_csv(ROOT/"data/processed/shell_sensitivity.csv")
    contrasts = pd.read_csv(ROOT/"data/processed/paired_contrasts.csv")
    audit = json.loads((ROOT/"results/interface_diagnostics/interface_audit.json").read_text())
    body = [r"""\documentclass[11pt,letterpaper]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,booktabs,tabularx,hyperref,xurl}
\hypersetup{hidelinks}
\begin{document}
\begin{center}
{\Large\bfseries S1 Text. Supporting methods and numerical results}\par
\medskip
Electrostatic--hydrophobic cooperation in the mesoscale assembly of synapsin and synaptic vesicles\par
\medskip
Shokoofeh Akbari, Yousef Jamali, and Oleg Shupliakov
\end{center}
\section*{S1. Design, estimands, and provenance}
All numerical tables refer to the 12 registered confirmatory trajectories:
four E/H conditions and seed blocks 9310002, 9310003 and 9310004.
Each run contains 200 SV and 1,807 Syn agents in a spherical compartment
of volume $0.1\,\mu\mathrm m^3$, evolved for 28 ms at a fixed integration
step of 0.0007 ms. Van der Waals interactions are enabled in all conditions.
Saved frames are separated by approximately 0.14 ms. Late-window
quantities are first averaged over frames 141--200 within a run (about
19.74--28 ms), then summarized across the three runs; agents and frames
are not independent replicates. Reported decimal values are rounded only
for presentation. Machine-readable full-precision tables are in S1 Data.

The association graph uses $\|x_i-x_j\|\leq a_i+a_j+2\,\mathrm{nm}$,
where $a_i$ is the outer radius. $L$ is the all-agent largest-component
fraction, $Q_{\rm Syn}$ is the non-isolated Syn fraction, and $L_{\rm SV}$
is the fraction of all SV in the all-agent largest component.
Non-isolated Syn includes Syn--Syn-only groups and is not a membrane-binding
measurement. The underlying data column ``synapsin\_recruited\_fraction''
retains its historical machine-readable name and means exactly $Q_{\rm Syn}$.
The primary outcome is $L$; the other two are complementary registered outcomes.

\section*{S2. Complete per-run registered outcomes}
"""]
    rows = []
    for _, r in endpoints.iterrows():
        rows.append([f"E{r.E}H{r.H}", str(r.seed),
            *[f"{r[c]:.6f}" for c in ("late_mean_largest_component_fraction",
                  "late_mean_synapsin_recruited_fraction", "late_mean_largest_component_vesicle_fraction")],
            f"{r.mean_drift_limited_fraction:.6f}"])
    body.append(table(["Arm", "Seed", "$L$", "$Q_{\\rm Syn}$", "$L_{\\rm SV}$", "Drift-limited"], rows, "llrrrr"))
    body.append(r"""The final column is the recorded whole-run mean drift-limited fraction,
not a fraction of Brownian increments clipped. The limiter affects deterministic
drift only. The runtime geometric screen is not proof of exact continuous-time
integration, equilibrium, or thermodynamic phase coexistence.
\clearpage
\section*{S3. Seed-block interaction and uncertainty}
The primary interaction is $I_s=L_{11,s}-L_{10,s}-L_{01,s}+L_{00,s}$.
It compares joint activation with an algebraic additive reference, not a fifth
simulated condition. All seed values are retained.
""")
    body.append(table(["Seed", "$I_s$"], [[k, f"{v:.9f}"] for k, v in audit["primary_interaction_by_seed"].items()], "lr"))
    c = contrasts[(contrasts.metric == "late_mean_largest_component_fraction") & contrasts.contrast.eq("E_H_interaction")].iloc[0]
    body.append(f"Mean primary interaction: {c['mean']:.6f}; between-seed SD: {c.sd:.6f}; approximate 95\\% paired Student-$t$ interval: [{c.ci95_low:.6f}, {c.ci95_high:.6f}].\n")
    body.append(r"""There are only three paired blocks (two degrees of freedom).
All eight sign-flip assignments produce a minimum possible two-sided
probability of 0.25 for nonzero effects; no high-resolution nonparametric
significance claim is made. Paired parametric intervals assume approximately
normal block differences, which cannot be established with three blocks.
Full conditional contrasts, exact sign-flip results, paired-$t$ estimates,
and separate within-outcome Holm adjustments are in S1 Data.

A supplementary ordinary-least-squares model includes intercept, E, H,
E-by-H and two seed indicators. Its design rank is six and residual degrees
of freedom are six. HC3 covariance addresses heteroskedasticity, not arbitrary
within-block dependence; these factorial intervals are secondary descriptive
checks rather than substitutes for additional independent trajectories.
An independently implemented matrix calculation reconciles all 15 released
paired-contrast means and interval endpoints to absolute tolerance $10^{-12}$.
\section*{S4. Association-shell sensitivity}
Only the graph threshold changes; the trajectories and forces do not.
Each entry below is a mean of three per-run late-window means.
""")
    perrun = shells[shells.frame_index >= 140].groupby(["E", "H", "seed", "analysis_shell_nm"])[["largest_component_fraction", "largest_component_vesicle_fraction"]].mean().reset_index()
    summary = perrun.groupby(["E", "H", "analysis_shell_nm"])[["largest_component_fraction", "largest_component_vesicle_fraction"]].mean().reset_index()
    rows = [[f"E{int(r.E)}H{int(r.H)}", f"{r.analysis_shell_nm:.0f}", f"{r.largest_component_fraction:.6f}", f"{r.largest_component_vesicle_fraction:.6f}"] for _, r in summary.iterrows()]
    body.append(table(["Arm", "Shell (nm)", "$L$", "$L_{\\rm SV}$"], rows, "lrrr"))
    body.append(r"""\clearpage
\section*{S5. Post hoc interface analysis}
Let $\mathcal S$ and $\mathcal V$ be the Syn and SV index sets. For each Syn,
\[
d_i^{\rm SV}=\sum_{j\in\mathcal V}\mathbf1\{\|x_i-x_j\|\leq a_i+a_j+2\,\mathrm{nm}\}.
\]
The two Syn fractions count $d_i^{\rm SV}\geq1$ and $d_i^{\rm SV}\geq2$,
respectively, divided by $N_{\rm Syn}$. The reciprocal SV coverage counts
SV with any Syn neighbor, divided by $N_{\rm SV}$. These are proximity
relations, not biochemical bonds or bridge lifetimes. This descriptive
analysis was added after the simulations and is outside the registered
testing family. No condition-specific graph threshold was tuned.

Complete per-run late-window fractions are given below.
""")
    rows = [[f"E{int(r.E)}H{int(r.H)}", str(int(r.seed)), f"{r.syn_near_sv_fraction:.6f}", f"{r.syn_multivesicle_neighbor_fraction:.6f}", f"{r.sv_near_syn_fraction:.6f}"] for _, r in interface.iterrows()]
    body.append(table(["Arm", "Seed", "Syn: $\\geq1$ SV", "Syn: $\\geq2$ SV", "SV: $\\geq1$ Syn"], rows, "llrrr"))
    body.append("Complete unordered pair counts (mean number of edges per snapshot within each run):\n")
    rows = [[f"E{int(r.E)}H{int(r.H)}", str(int(r.seed)), f"{r.edges_SS:.3f}", f"{r.edges_SV:.3f}", f"{r.edges_VV:.3f}"] for _, r in interface.iterrows()]
    body.append(table(["Arm", "Seed", "Syn--Syn", "Syn--SV", "SV--SV"], rows, "llrrr"))
    body.append(r"""The spatial-tree graph was compared with all-pairs enumeration on toy
and random geometries and on the first and last saved frame of every run
(24 trajectory frames). Boundary inclusion, translation, permutation and
no-edge tests passed. Reconstructed $L$, $Q_{\rm Syn}$ and $L_{\rm SV}$
matched all 12 released per-run tables to absolute tolerance $10^{-12}$.
Input HDF5 checksums matched the availability manifest.
\clearpage
\section*{S6. Post hoc mobility and the experimental scale reference}
Syn positions are projected onto $xy$, converted from nm to $\mu$m,
and fitted using per-agent time-averaged MSD with a free intercept over
native lags 2.10--3.92 ms (14 lag points). The reported apparent coefficient
is slope/4 in $\mu\mathrm m^2/\mathrm s$. Negative fitted slopes are retained.
The whole 200-frame trajectory is the principal scale diagnostic; the final
60-frame window is a sensitivity check, not a result-selected replacement.
""")
    rows = []
    for _, r in mobility[mobility.window == "full_28ms"].sort_values(["arm", "seed"]).iterrows():
        late = mobility[(mobility.arm == r.arm) & (mobility.seed == r.seed) & mobility.window.eq("late_60frames")].iloc[0]
        rows.append([r.arm, str(r.seed), f"{r.D_mean_um2_s:.6f}", f"{late.D_mean_um2_s:.6f}", f"{r.negative_slope_fraction:.4f}", f"{late.negative_slope_fraction:.4f}"])
    body.append(table(["Arm", "Seed", "$D$: full", "$D$: late", "Negative: full", "Negative: late"], rows, "llrrrr"))
    body.append(r"""Stationary and constant-velocity analytic tests passed. The mean
per-agent estimate matches an independently fitted ensemble-MSD slope to
$10^{-12}$. IDs, species, frame counts and trajectory hashes are verified.

Hoffmann et al. (Nature Communications 2023, DOI
\href{https://doi.org/10.1038/s41467-023-42372-6}{10.1038/s41467-023-42372-6})
provide a Fig. 2h short-time mobility scale: source-data means are
0.0189192370 with SV and 0.0253090377 without SV, in
$\mu\mathrm m^2/\mathrm s$. Source data:
\url{https://doi.org/10.6084/m9.figshare.24179445.v1}.
No experimental uncertainty is reconstructed. Their experiment selects
condensate-core tracks at a 1-ms cadence; this model selects all Syn at
a 0.14-ms saved cadence without localization error or exposure blur.
Every model arm contains SV. Furthermore, the effective friction was
previously literature-informed. The comparison is consequently a scale
diagnostic, not matched perturbation validation or independent parameter
calibration. The absolute model mobility remains higher.

\section*{S7. Reproduction and data interpretation}
S1 Data contains full-precision numerical inputs for the four figures,
registered tables and post hoc diagnostics, plus an inventory with SHA-256
checksums. The article-specific code package contains standalone extraction,
statistics, rendering and supporting-table scripts. Published input
trajectories are enumerated in its availability manifest; original raw
simulation outputs are preserved separately. No trajectory is discarded
or modified by the post hoc diagnostics.

These results describe finite-horizon structural organization in the
implemented model. They do not establish stationary coexistence, molecular
bond kinetics, disease causality or a unique pair-class mechanism. Whole-kernel
switches alter several pair types; proximity configurations identify useful
structural readouts rather than isolated molecular causes.
\end{document}
""")
    atomic_text(a.manuscript/"S1_Text.tex", "\n".join(body))
    files = list((ROOT/"data/processed").glob("*.csv"))
    files += list((ROOT/"results/interface_diagnostics").glob("*.csv"))
    files += list((ROOT/"results/mobility_diagnostic").glob("*.csv"))
    files += list((ROOT/"results/statistical_audit").glob("*.csv"))
    files += [ROOT/"data/external/hoffmann_mobility_reference.json", ROOT/"DATA_LICENSE.md"]
    files += [ROOT/f"results/{d}/{f}" for d, f in [("interface_diagnostics", "interface_audit.json"),
              ("mobility_diagnostic", "audit.json"), ("statistical_audit", "audit.json")]]
    inventory = [{"file": str(f.relative_to(ROOT)).replace("\\", "/"),
                  "bytes": f.stat().st_size, "sha256": hashlib.sha256(f.read_bytes()).hexdigest()}
                 for f in sorted(files)]
    temporary = a.output/".S1_Data.editing.zip"
    with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for f, record in zip(sorted(files), inventory):
            z.write(f, record["file"])
        z.writestr("INVENTORY.json", json.dumps(inventory, indent=2))
        z.writestr("README.txt", "S1 Data: numerical inputs for the 12-run force study.\n"
            "Frames/agents are nested, not independent replicates.\n"
            "Historical column synapsin_recruited_fraction means non-isolated Syn fraction.\n"
            "Interface and mobility diagnostics are post hoc; see audit JSON files.\n"
            "Author-generated data: CC BY 4.0. Experimental values retain original terms.\n")
    os.replace(temporary, a.output/"S1_Data.zip")
    with zipfile.ZipFile(a.output/"S1_Data.zip") as z:
        if z.testzip() is not None:
            raise ValueError("ZIP integrity failure")
        for record in inventory:
            if hashlib.sha256(z.read(record["file"])).hexdigest() != record["sha256"]:
                raise ValueError("ZIP content checksum failure")
    print(f"S1 Text generated; {len(files)} numerical inputs verified in S1_Data.zip")


if __name__ == "__main__":
    main()
