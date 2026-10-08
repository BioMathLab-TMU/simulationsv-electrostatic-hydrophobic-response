# Reproduce the article's additional diagnostics and figures

Run from the standalone repository after installing `requirements.txt`.
The complete released input trajectories are at
`data/raw/<run_id>/trajectory.h5`. Each extractor verifies its inputs against
`data/availability.json`.

```bash
python analysis/interface_diagnostics.py --output results/interface_diagnostics
python analysis/mobility_diagnostic.py --output results/mobility_diagnostic
python analysis/statistical_audit.py
python -m unittest discover -s tests -v
python scripts/render_manuscript_figures.py --output figures
python scripts/prepare_supporting_data.py --manuscript build/supporting --output build/submission
```

The final command creates a standalone supporting-methods source file and
S1 Data ZIP. No manuscript checkout or private writing workspace is required.
The analyses do not modify any simulation configuration, trajectory or
registered endpoint.

## Observation definitions

At shell=2 nm, an edge satisfies Euclidean center distance <= outer radius sum
plus shell. Unordered pairs are counted once. Syn adjacent to at least two SV
is a geometric bridge candidate, not a molecular bond or a lifetime. The
registered historical data label `synapsin_recruited_fraction` denotes the
fraction of Syn with any neighbor, including Syn-only components.

Mobility uses xy projection, nm-to-um conversion, seconds as fit time, a free
intercept and slope/4 over 2.10--3.92 ms. Full and late windows are both retained,
as are negative individual fitted slopes. The experimental reference is an
unmatched scale comparison, not an independent calibration or prediction test.

All additional interface/mobility readouts are post hoc. Their purpose is to
describe and interpret the existing factorial results, not to enlarge the
registered confirmatory testing family. Replication remains n=3 independent
seed blocks per condition; frame and agent counts do not enlarge that n.

## Supporting data and access

The generated `S1_Data.zip` contains numerical figure/table inputs and a SHA-256
inventory, not large raw HDF5 trajectories. Those trajectories are included
directly in this repository and identified by `data/availability.json`.
Editors and reviewers must have access at submission; a private URL alone is
insufficient. These analysis commands do not change repository visibility or
submit a manuscript.

Original software is MIT-licensed; original generated data/tables/figures are
CC BY 4.0. Third-party source data retain the provider's terms.
