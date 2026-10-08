# Released data dictionary

## Independent units and identities

One trajectory/run is one independent replicate. E/H conditions share three
seed-block identifiers (9310002, 9310003, 9310004). Frames, graph components
and particles are nested observations and do not increase the replicate count.
Species codes: 1 = Syn, 2 = SV. IDs are stable through saved frames.

## HDF5 input

Each run folder contains `trajectory.h5`, `Dashboard.input.json` and
`numerical_diagnostics.json`. File hashes and config hashes are in
`data/availability.json` and the design manifest. Frame groups are
`step_200` through `step_40000` at intervals of 200 integration steps.
Each contains 2,007 agents and an explicit physical timestamp.

| Field | Meaning / units |
| --- | --- |
| positions | Shape (2007, 3), center coordinates, nm |
| radii, core_radii | Outer effective and membrane/core radius, nm |
| ids, types | Stable particle identity and species code |
| physical_time_ms | Time of the recorded state, ms |
| simulation_step | Physical integration-step index |
| delta_t_ms, frame_interval_ms | Integration step and scheduled saved cadence, ms |
| *_force_vector | Per-agent force components in internal force units |
| total_force_magnitude | Magnitude in internal force units |
| drift_scale | Multiplier applied to deterministic drift, not Gaussian noise |
| previous_step_overlap_projection_* | Diagnostics of the preceding correction |
| frame_write_complete | Saved-state completion indicator |

One internal force unit is 1.380649e-14 N. See `docs/UNITS.md` for all
conversion conventions. Force values and positions refer to the recorded
state; previous-step projection diagnostics must not be relabelled as
simultaneous force measurements. `run_metadata`, `run_completion` and
`minimizer_diagnostics` retain the selected numerical metadata.

## Registered outcomes

At a surface-gap shell, an unordered pair is adjacent if center distance is
at most the sum of outer radii plus shell. The primary shell is 2 nm.

- `largest_component_fraction`: size of the all-agent largest component / 2007.
- `synapsin_recruited_fraction`: historical label for **non-isolated Syn**
  / all Syn. This includes Syn-only groups, not only Syn-SV association.
- `largest_component_vesicle_fraction`: SV in the all-agent largest
  component / all 200 SV, not a separately optimized SV-only component.
- Prefix `late_mean_`: last-60-frame mean within a run, before summarizing
  across three seeds. No framewise pooling across conditions is used.

Complete time series, shell sensitivity and run-level estimates are in
`data/processed/`. Counts need not be integers after averaging within a run.
Mean/SD and uncertainty calculations must use per-run values.

## Post hoc outputs

`results/interface_diagnostics/` contains per-frame and per-run SS/SV/VV
unordered edge counts, Syn with one or multiple SV neighbors and reciprocal
SV coverage. These are geometric proximity relations, not molecular bonds.

`results/mobility_diagnostic/` contains xy-TAMSD values, individual-agent
fitted apparent coefficients and run summaries for full/late windows. Position
is converted to micrometers and fit time to seconds; coefficient units are
micrometer squared per second. Negative fitted slopes are retained.

`results/statistical_audit/` contains independently reconstructed paired
contrasts and approximate OLS/HC3 outputs. With n=3, a two-sided exact
sign-flip probability cannot be below 0.25 before multiplicity adjustment.
The HC3 calculation is not a cluster-robust small-block procedure.

## Provenance and exclusions

The release preserves all values in the included HDF5 datasets exactly.
Original machine paths, RNG checkpoint state and unrelated diagnostics are
not needed for these reported analyses and are excluded from the public
subset. Original full simulation outputs are not deleted or overwritten.
No personally written project reports or other-article datasets are released.
