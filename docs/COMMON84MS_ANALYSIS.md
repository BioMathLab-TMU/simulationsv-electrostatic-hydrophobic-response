# Common 84-ms analysis

The four electrostatic/hydrophobic conditions are observed for 84 ms in
three seed blocks. Physical settings, particle identities and numerical
settings are unchanged from the original 28-ms trajectories. Checkpoint
joins were verified against positions, radii and RNG state, then duplicate
join frames were counted once. There are 600 unique frames per trajectory.

The main structural window is the last 30%: frames 421--600, or
58.94--84 ms. The historical 28-ms window remains frames 141--200,
19.74--28 ms. These are repeated observations of the same trajectories,
with three independent seed blocks. Frames and agents are not replicates.

The graph threshold remains distance <= sum of outer radii + 2 nm.
Sensitivity shells are 0, 4 and 6 nm. Each unordered pair is counted once.
The largest component is selected by all-agent size; its SV fraction
uses the number of all SVs as denominator. Non-isolated Syn includes
Syn--Syn groups. Contact categories resolve proximity, not molecular bonds.

The apparent projected mobility uses per-agent xy TAMSD and free-intercept
OLS over 14 native lags, 2.10--3.92 ms. Slopes are divided by four.
The full 600-frame and final 180-frame estimates remain separate.

All raw segments passed their recorded numerical criteria. Saved geometry
is float32 while the solver uses double precision. One reconstructed
residual exceeded 0.005 nm by 1.52e-6 nm; its value lies within the
per-coordinate half-ULP rounding bound, and the recorded solver residual
was below 0.005 nm. This does not change the solver tolerance.

The complete earlier 2,400-frame structural/contact results were recovered
within 1e-12. Graph calculations on selected frames were also compared
with all-pairs enumeration. Running analysis/common84ms.py --verify
reproduced the current structural, contrast and mobility tables from the
public analysis-input HDF5 files. Original raw files remain unchanged.

The later data retain the largest assembly under joint activation and
show slower late growth in that condition. The separate interactions
continue to build networks, reducing the earlier interaction contrast.
At 84 ms, the mean all-agent contrast is 0.110307 with mixed seed signs;
the mean SV-participation contrast is -0.050037. These are finite-time
structural comparisons. The repository does not certify equilibrium.

Current data: data/common84ms_availability.json and
data/processed/common84ms/. Original 28-ms data remain at their previous
paths. Both versions are retained by separate manuscript tags.
