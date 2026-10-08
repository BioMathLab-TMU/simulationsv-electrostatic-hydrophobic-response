# SimulationSV: electrostatic-hydrophobic response

Article-specific reproducibility materials from the Biomathematics Laboratory,
Tarbiat Modares University. This repository supports the study
**Electrostatic-hydrophobic cooperation in the mesoscale assembly of synapsin
and synaptic vesicles** by Shokoofeh Akbari, Yousef Jamali, and Oleg Shupliakov.

Repository: https://github.com/BioMathLab-TMU/simulationsv-electrostatic-hydrophobic-response

This is a code/data release for a research manuscript, not a claim of journal
acceptance. Results concern **finite-time structural organization in the
implemented model**, not established equilibrium or thermodynamic LLPS.

## Study and released data

- Fixed 2 x 2 electrostatic (E) / hydrophobic (H) design; van der Waals
  interactions remain enabled in all four conditions.
- Three independent seed blocks per condition: 12 trajectories in total.
- 200 vesicles and 1,807 effective synapsin particles; 30 micromolar nominal
  synapsin concentration in a 0.1 cubic-micrometer spherical compartment.
- Fixed timestep 0.0007 ms, 28-ms horizon, 200 saved frames per trajectory.
- Primary 2-nm association graph; observation-shell checks at 0, 4 and 6 nm.
- Interface and projected-mobility diagnostics are explicitly post hoc.

All 12 released HDF5 files are in `data/raw/`, approximately 606 MB in total.
They contain positions, stable IDs, species, radii, times, force vectors and
selected numerical diagnostics. They are not lossy resamplings: included
datasets preserve the original values exactly. Machine provenance, RNG
checkpoint state and unrelated channels are excluded; original simulation
outputs remain preserved separately. Configurations and data SHA-256 hashes
are recorded in the manifests.

## Where to start

| Location | Content |
| --- | --- |
| Root C++ source, CMakeLists.txt, third_party/ | Compilable simulator and dependency notices |
| configs/manifest.json, configs/runs/ | The exact 12-run design and checksummed inputs |
| data/availability.json, data/raw/ | Completed input trajectories and availability metadata |
| data/processed/ | Registered outcomes, per-frame values, paired contrasts, shell checks |
| results/interface_diagnostics/ | Complete post hoc interface values and audit |
| results/mobility_diagnostic/ | Per-run/per-agent mobility estimates, MSD inputs and audit |
| results/statistical_audit/ | Independent paired-contrast and HC3 calculations |
| figures/ | Registered checks and the four manuscript figure sets |
| docs/ | Equations, units, parameter provenance, data dictionary and reproduction guide |
| provenance/ | Source/data checksums, release inventory and tested environment |
| tests/ | Unit, geometry, data-integrity and statistical-contract tests |

No personal project reports, other-study data, manuscript drafts, cover letters,
private writing history or machine-specific build folders are part of this release.

## Reproduce the analysis without running new simulations

Python 3.12 is the tested interpreter. Linux and Windows/WSL are supported.
Use an isolated environment and install the declared dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/reproduce.py --input data/raw --output results/recomputed --workers 1 --require-complete
python analysis/interface_diagnostics.py --output results/interface_diagnostics
python analysis/mobility_diagnostic.py --output results/mobility_diagnostic
python analysis/statistical_audit.py
python scripts/render_manuscript_figures.py --output figures
python scripts/prepare_supporting_data.py --manuscript build/supporting --output build/submission
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1`.
The last command generates supporting-methods LaTeX and S1 Data numerical
inputs, not a new manuscript. It does not require any private project checkout.
Figure rendering prefers Arial; a declared available fallback may be used on
Linux. Released journal figure files use embedded Arial and 600-dpi LZW TIFF.
See `docs/ARTICLE_DIAGNOSTICS.md` for estimands and comparison limitations.

The independent unit is a complete run/seed block, **not an agent or frame**.
There are only three blocks per arm. Proximity does not establish biochemical
bonds or contact lifetimes. The experimental mobility comparison is an unmatched
scale diagnostic; the absolute model mobility is higher.

## Build and rerun the physical model (optional)

A C++23 compiler, CMake, OpenMP, Boost and HDF5 development libraries are needed.
For Ubuntu/WSL:

```bash
sudo apt-get install build-essential cmake libboost-system-dev libboost-filesystem-dev libhdf5-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
python scripts/run_study.py --executable build/SimulationSV
```

This optional command launches new simulations. It is not needed to reproduce
the published analysis. Runs execute sequentially; the default runtime policy
reserves two detected physical cores. Consult `docs/METHODS.md` before interpreting
the implemented time scale, drift limiter and geometric correction.

## External source data

Hoffmann et al. (2023), DOI 10.1038/s41467-023-42372-6, provide the experimental
mobility reference. The original source data are at
https://doi.org/10.6084/m9.figshare.24179445.v1 .
This repository includes extracted scale-reference values and provenance,
not the publisher PDF or full third-party workbook. Original provider terms apply.

## Authors, correspondence and funding

Shokoofeh Akbari and Yousef Jamali: Department of Applied Mathematics, School of
Mathematical Sciences, Tarbiat Modares University, Tehran, Iran.
Oleg Shupliakov: Department of Neuroscience, Karolinska Institutet, Stockholm, Sweden.
Correspondence: Yousef Jamali, y.jamali@modares.ac.ir;
ORCID https://orcid.org/0000-0002-5709-3821 .

This work is based upon research funded by Iran National Science Foundation
(INSF) under project No. 4029556. The funder had no role in study design, data
collection and analysis, decision to publish, or preparation of the manuscript.

## Licensing and citation

Original author-generated code and software documentation: MIT (`LICENSE`).
Original simulation data, tables and research figures: CC BY 4.0
(`DATA_LICENSE.md`). Third-party components retain their original notices.
Use `CITATION.cff` and identify the exact release tag or commit used. No article
DOI is fabricated; publication details will be added when available.
