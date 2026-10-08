# Third-party dependencies

The bundled `highfive/` headers retain their original copyright headers and
are distributed under Boost Software License 1.0. See
`licenses/Boost-1.0.txt` and the upstream project:
https://github.com/BlueBrain/HighFive

The bundled single-header `nlohmann/json.hpp` is version 3.11.3, copyright
Niels Lohmann and contributors, licensed under MIT. See
`licenses/nlohmann-json-MIT.txt` and https://github.com/nlohmann/json .

`Srock2Coefficients.h` contains coefficient code with the
StochasticDiffEq.jl MIT attribution preserved in
`licenses/StochasticDiffEq-MIT.txt`. This experimental integrator is compiled
but not selected by any registered run; every run uses Euler–Maruyama.

The Boost, OpenMP and HDF5 system libraries are not vendored into this
repository. Their own licenses apply to installations on the reader's machine.
