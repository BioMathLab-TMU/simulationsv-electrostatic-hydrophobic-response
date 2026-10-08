# Parameter provenance and limits

The exact numerical inputs are the 12 JSON files under `configs/runs/`; their
SHA-256 hashes are registered in `configs/manifest.json`. The generated
`docs/parameters.csv` lists the common and arm-varying input values with units
and their evidentiary status. A source supporting the *scale* or *sign* of a
coarse-grained parameter is not misrepresented as measuring the implemented
effective parameter directly.

| Parameter choice | Evidentiary status and primary source |
| --- | --- |
| 0.1 µm³ spherical bouton; 200 SVs | Chosen within measured order of magnitude from [quantitative 3-D EM of hippocampal boutons](https://pmc.ncbi.nlm.nih.gov/articles/PMC6573206/). Sphere is a modeling geometry, not a claim that every bouton is spherical. |
| SV core mean 20.45 nm, SD 2.5 nm | Half the reported 40.9 ± 5.0 nm diameter in [cryo-electron tomography](https://pmc.ncbi.nlm.nih.gov/articles/PMC5815350/). The 15–30 nm truncation is a defensible modeling bound consistent with [SV morphometry](https://pmc.ncbi.nlm.nih.gov/articles/PMC2849067/). |
| SV effective corona +4.8 nm | Outer-versus-bilayer radius difference represented from [SV structure/intervesicle-spacing analysis](https://pub.dzne.de/record/273920/files/DZNE-2024-01394.pdf?subformat=pdfa); a uniform spherical shell is an abstraction. |
| Syn effective radius 3.972 nm | Fixed effective modeling radius, not the experimentally measured radius of gyration. [Synapsin structural measurements](https://pmc.ncbi.nlm.nih.gov/articles/PMC11628805/) constrain the scale, not exact equality. |
| 30 µM Syn, 1807 agents | Study-selected concentration, with count derived from volume and Avogadro's constant; not an experimentally measured bouton count. |
| Temperature 310.15 K; water εr 74.2 | Physiological temperature and a temperature-adjusted water prior; [NIST water-permittivity measurements](https://www.nist.gov/publications/static-dielectric-constant-liquid-water-between-274-and-418-k-near-saturated-vapor-0) support its order of magnitude. |
| Effective viscosity 0.35 Pa·s | Coarse-grained mobility parameter. It is **not** a direct free-water viscosity; associated kinetic time scales remain model-dependent. |
| I=0.15 M; λD≈0.779 nm | Study-selected physiological-screening condition; Debye length derived from SI constants and εr rather than fitted to the outputs. |
| Syn +16 mV, SV −22.5 mV | Effective electrostatic parameters. [Synapsin positive-charge/multivalency evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC8424992/) and [negative SV zeta-potential measurements](https://pmc.ncbi.nlm.nih.gov/articles/PMC2852195/) support sign/scale; neither is a direct joint measurement of these exact potentials in the simulated geometry. |
| Hydrophobic SS 1.5 kBT, SV 1.0 kBT; decay 0.5 nm | Fixed model-prior contact strengths and short-range decay. [Native-like SV-binding experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC2115549/) support combined hydrophobic/electrostatic binding but do not uniquely identify these pair energies. The 0.5-nm decay is a coarse-grained short-range choice, not a directly measured Syn–SV decay constant. |
| Hamaker A=5×10⁻²¹ J, minimum separation 0.3 nm | Fixed nonspecific background from [lipid-bilayer force measurements](https://www.jstage.jst.go.jp/article/jos/67/11/67_ess18088/_html/-char/en); kept ON in all arms. Different material pairs are approximated by this common effective value. |
| dt=0.0007 ms, 128 projection ceiling, 0.5-nm core allowance, R/4 drift cap | Numerical design decisions of this discrete model, not direct measurements from biology. Recorded QC assesses admissibility under the protocol, not global solver convergence. |
| 2-nm primary edge shell; 0/4/6-nm sensitivity | Prespecified graph-analysis definitions. They do not impose a 6-nm physical-force cutoff. |

The four arm controls change **only** E/H switches; all 12 configurations can
be compared directly against the registered manifest. Van der Waals forces
are still calculated, but this article does not allocate an independent
ablation or infer their causal role from this 2×2 design.
