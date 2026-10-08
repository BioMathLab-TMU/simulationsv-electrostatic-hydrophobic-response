# Processed article results

All files in this directory come from `scripts/reproduce.py` applied to the
12 registered, included runs. The `run_inventory.csv` is the inclusion/QC
decision for every planned run. `run_endpoints.csv` contains the independent
seed-level values. `frame_metrics.csv` and `shell_sensitivity.csv` support
re-plotting and threshold checks; frames are not statistical replicates.

`condition_summary.csv` reports arm means, SDs and 95% t intervals.
`paired_contrasts.csv` contains the five prespecified contrasts on each
endpoint. The sign-flip p-values are descriptive at n=3 (minimum 0.25);
Holm-adjusted values are given. None of these p-values licenses a blanket
claim about all possible parameter settings, equilibration, or LLPS.

The figures in `../../figures/` use these tables. The per-frame outputs and
figures can be regenerated from `../raw/` with the command in the root README.
Output paths, cache metadata and machine-local identifiers are excluded from
these released tables.
