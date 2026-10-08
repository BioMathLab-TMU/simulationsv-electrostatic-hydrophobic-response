# Study data

Each `raw/<run_id>/` directory contains one immutable, numerically complete
trajectory subset (`trajectory.h5`), the exact registered input configuration
(`Dashboard.input.json`), and numerical diagnostics. All 12 registered
condition/seed combinations are listed in `availability.json`, along with
SHA-256 hashes and the number of saved frames.

The HDF5 export retains all saved positions, types, IDs, outer/core radii,
physical times, force components and totals, drift-scale and projection
diagnostics, initialization data, model settings, and run-completion data.
Every included dataset was compared element-by-element (including dtype and
shape) to the completed original before publication. There is no interpolation,
downsampling, or precision conversion of an included field.

The export intentionally excludes the original machine-provenance group, RNG
checkpoint-state bytes, and unrelated scalar summaries of pair-force
distributions. Those original full HDF5s remain in the private working
archive; this public repository contains exactly what is needed for the
current article's planned analysis. The excluded data are not used by
`scripts/reproduce.py`.

The HDF5 files are stored directly in Git, so a normal clone obtains the full
trajectory data without a separate large-file transfer. The released files are
immutable; do not edit the HDF5 trajectories in place.
`scripts/export_completed_data.py` documents how the reduced, value-verified
release was generated from the completed originals.
