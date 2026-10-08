# Registered analysis protocol

## Design and independence

The confirmatory design consists of three independent seeds (9310002–9310004),
each used in the four conditions E1H1, E1H0, E0H1, and E0H0. E and H denote
electrostatic and hydrophobic force switches. Using the same seed across arms
defines a paired block; trajectories need not remain identical after the force
laws diverge. The preliminary seed 9310001 is not part of confirmatory inference.

Every arm uses 200 vesicles and 1807 synapsin agents in a 0.1 µm³ spherical
domain, the same radius law, the same timestep, the same 28-ms horizon, and the
same numerical controls. Van der Waals interactions remain enabled.

## Geometric observables

At each saved frame, nodes are particles and an undirected edge is present when
the direct Euclidean surface gap satisfies

`distance(x_i,x_j) - R_outer,i - R_outer,j <= shell_nm`.

The primary shell is 2 nm. The prespecified sensitivity shells are 0, 4, and
6 nm. These are **analysis definitions**, not new force cutoffs. No periodic
minimum-image operation is used. Components are found over all particles.
The largest-component fraction is its node count divided by the total number
of particles. Synapsin recruitment is the fraction of synapsin nodes in
components of size at least two. The largest-component vesicle fraction uses
all vesicles as its denominator, not all particles.

The primary endpoint is the mean all-particle largest-component fraction over
the final 30% of the 200 saved frames (60 frames, 19.74–28 ms). Secondary
endpoints are the analogous synapsin-recruitment and largest-component-vesicle
fractions. Equal saved-frame spacing makes this an arithmetic time-window
mean. The late window is not called equilibrium.

## Paired contrasts

For each seed and endpoint, use its four arm values to form

* E effect at H=0: `Y10 - Y00`;
* E effect at H=1: `Y11 - Y01`;
* H effect at E=0: `Y01 - Y00`;
* H effect at E=1: `Y11 - Y10`;
* interaction: `Y11 - Y10 - Y01 + Y00`.

The last quantity is an additive-scale factorial interaction, not an assertion
of biological synergy independent of the chosen observable and model.
Report all seed-level values, means, SD, and paired mean-difference 95% t
intervals with two degrees of freedom. These intervals are approximate and
unstable at n=3; do not hide their width or clamp them to a fraction's bounds.
Exact two-sided sign-flip p-values are included only with their symmetry
assumption and coarse resolution stated. At n=3 their smallest attainable
value is 0.25. Holm-adjusted p-values use the five contrasts within each
endpoint family; secondary endpoints and shell sensitivity remain secondary.
No frames or particles are promoted to independent replicates.

## Inclusion and limitations

All registered runs are accounted for. A complete, numerically admissible run
with weak, absent, or opposite structural response remains a valid observation.
Only prespecified numerical/integrity failures are excluded; reasons are saved.
Missing data are reported, never filled with desired outcomes. Final inference
requires all 12 registered runs, correct seeds/configurations, the full saved
schedule, finite states, and runtime-QC acceptance.

Numerical QC checks completion and recorded admissibility. It is not proof of
timestep convergence, equilibrium, contact-event resolution, or LLPS. The
drift cap and projection define the implemented discrete update and must be
reported, including how often they activate. The main conclusions must remain
within finite-time model organization and the tested parameter settings.
