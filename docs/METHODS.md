# Model and numerical methods

Synapsin and synaptic vesicles are effective spherical agents with species-
specific radii and mobilities. The implementation uses position-dependent
electrostatic, hydrophobic, and van der Waals pair forces plus Brownian noise.
The study is a coarse-grained organization model, not a molecular-resolution
protein-conformation or complete exocytosis/endocytosis model.

## Geometry

Vesicle membrane/core radii are drawn from a normal law with mean 20.45 nm and
SD 2.5 nm, truncated to [15,30] nm. An effective 4.8-nm protein corona is added:
`R_outer = R_core + 4.8 nm`. Synapsin's effective radius is 3.972 nm. Corona
interpenetration is allowed; bilayer crossing is not treated as physical.
The numerical core allowance is 0.5 nm, with a 0.005-nm projection residual
tolerance. The permitted outer-radius overlap is consequently 0.5 nm for SS,
5.3 nm for SV, and 10.1 nm for VV, plus the stated numerical tolerance.

The spherical domain volume is 0.1 µm³ = 10⁸ nm³. Its radius is
`(3V/(4π))^(1/3) = 287.941191 nm`; its containing cube has side
575.882382 nm. A linked-cell side of 71.985298 nm exceeds the maximum outer
vesicle diameter 69.6 nm. The force search uses the fixed two-cell neighborhood.
All candidate pairs found in that neighborhood are evaluated without an extra
radial cutoff. This is **not an infinite-range/all-pairs force summation**:
the neighborhood remains the spatial truncation of the implemented model.
The 6-nm analysis shell has no separate role as a force switch.

The boundary is closed, reflecting, and nonperiodic. It describes confinement
within an idealized domain, not exchange through an active zone. The C++ source
is authoritative for reflection and position-correction details.

## Force kernels and units

Electrostatics use a linearized-Poisson–Boltzmann/Yukawa sphere interaction with
surface potentials +16 mV (Syn) and −22.5 mV (SV), at 0.15 M ionic strength.
The selected Debye length is 0.778856376 nm. Hydrophobic attraction uses an
exponential contact PMF with decay length 0.5 nm and contact energies 1.5 kBT
(SS), 1.0 kBT (SV), and zero (VV). Van der Waals attraction uses a nonretarded
Hamaker/Derjaguin kernel, A=5×10⁻²¹ J, and a 0.3-nm minimum separation. Its SV
surface reference is the membrane/core, not the outer protein corona.
Exact force signs, regularization, and pair-radius handling are in
`Interaction.cpp` and `docs/parameters.csv`; simplified prose is not a replacement
for those kernels.

The internal bases are nm, ms, and a chosen mass scale of 1.380649×10⁻¹¹ kg,
as defined in `Converter.cpp`. Thus 1 internal energy unit is
1.380649×10⁻²³ J and 1 internal force unit is 1.380649×10⁻¹⁴ N.
The internal viscosity unit is 13.80649 Pa·s. Report SI parameters
using the explicit converter rather than interpreting internal numeric fields
as SI values. Metadata names with `_nm`, `_ms`, `_kbt`, or `_mV` retain the
explicit units indicated by those names. See `docs/UNITS.md` for the exact
numerical conversion.

## Discrete update

For a particle with drag γ, the Brownian diffusion relation is D=kBT/γ.
An Euler–Maruyama proposal adds deterministic force drift and a Gaussian
increment with covariance `2D dt I`. The deterministic drift is capped at one
quarter of the outer radius; Brownian noise is not scaled by this cap. Boundary
handling and mobility-weighted overlap projection follow the proposal.
Projection has a maximum of 128 iterations and exits early when converged.
These operations form the actual numerical model; a continuous SDE alone
would not completely specify the implementation.

The fixed timestep is 0.0007 ms. There are 40,000 movement steps, with snapshots
at steps 200,400,…,40000 (200 frames, separated by 0.14 ms). A short geometric
initialization minimizer precedes dynamics and is recorded separately. All
force-toggle comparisons keep numerical controls identical.

Runtime diagnostics record nonfinite proposals, boundary rejection, projection
residuals/nonconvergence, drift-limit activation, and extreme Brownian/proposal
increments. `analysis/quality_control.py` states the exact acceptance rules.
Saved-frame spatial checks additionally inspect the radius-dependent overlap
allowance. Successful QC means acceptance under these rules, not universal
numerical accuracy or absence of bias.
