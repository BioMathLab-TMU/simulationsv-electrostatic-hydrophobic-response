# Internal units and SI recovery

The converters are defined by `Lu=10^9`, `tu=10^3`, `Mu=1.380649e-11`, and
`chu=1.602e-19`. These are the actual code constants, not a conventional amu
system. Taking L₀=10⁻⁹ m, T₀=10⁻³ s, M₀=Mu kg, and Q₀=chu C gives:

| Quantity | One internal unit in SI | Recovery |
| --- | --- | --- |
| Position/radius | 10⁻⁹ m | internal × 10⁻⁹ |
| Time | 10⁻³ s | internal × 10⁻³ |
| Energy | 1.380649×10⁻²³ J | internal × kB(SI) |
| Force | 1.380649×10⁻¹⁴ N | internal × kB(SI)/L₀ |
| Viscosity | 13.80649 Pa·s | internal × M₀/(L₀ T₀) |
| Diffusion | 10⁻¹⁵ m²/s | internal nm²/ms × 10⁻¹⁵ |
| Electric charge | 1.602×10⁻¹⁹ C | internal × Q₀ |
| Permittivity | Q₀²/(energy₀ L₀) F/m | internal × this factor |

For example, eta=0.025350396806139718 corresponds to 0.35 Pa·s and
hamaker=362.148525801996 corresponds to 5×10⁻²¹ J. The converted Boltzmann
constant is 1 energy-unit/K; `kBT` therefore equals T in internal energy units.
A parameter with the suffix `_kbt` is already an energy divided by kBT, not
an unscaled internal energy. Volume in µm³ multiplies by 10⁹ to give nm³.
The 0.1-µm³ volume and 30-µM concentration give
`round(30e-6 × 0.1e-15 × N_A) = 1807` synapsin agents.

The Debye-length helper uses the exact SI elementary charge 1.602176634e-19 C
and vacuum permittivity 8.8541878128e-12 F/m. The converter retains its historical
rounded charge 1.602e-19 C and vacuum permittivity 8.854e-12 F/m. Both conventions
are recorded rather than silently changed. This reproducibility package does
not change the production equations or their conversion constants.

`D_syn`, `D_ves`, `drag_syn`, and `drag_ves` in the input are dimensionless
multipliers, not measured SI diffusion coefficients. Effective mobility is
computed from viscosity and particle radius in the C++ model. The selected
effective viscosity is a coarse-grained modeling choice and makes physical
kinetic interpretation conditional on its calibration.
