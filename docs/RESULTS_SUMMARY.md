# Registered force-response results and interpretation limits

This is a concise guide to the machine-readable result tables, not a manuscript
or a substitute for journal peer review. All 12 registered confirmatory runs
completed; numerical and saved-state inclusion checks passed. Each arm has
three independent seeds. The preliminary pilot seed is excluded.

The primary endpoint is the mean all-particle largest-component fraction over
the final 60 saved frames at a 2-nm surface-gap shell. Arm means (SD) are:

| Electrostatic | Hydrophobic | Mean (SD) |
| --- | --- | --- |
| off | off | 0.02376 (0.00174) |
| off | on | 0.09898 (0.00710) |
| on | off | 0.10529 (0.02091) |
| on | on | 0.55277 (0.08024) |

The prespecified additive-scale interaction is 0.37226, with a paired-seed
95% t interval [0.20676, 0.53776]. Its exact two-sided sign-flip p-value is
0.25 with n=3 (Holm-adjusted p=1.0 across the five primary contrasts). This
small-n result is a reproducible *model pattern across the sampled seeds*,
not strong nonparametric evidence at a 0.05 threshold. The interval relies on
a different, small-sample parametric approximation and should not be used to
hide the low resolution of the exact test.

The sign of E×H is endpoint-dependent: the interaction for synapsin recruitment
is −0.02893, whereas the one for the fraction of vesicles in the largest
component is +0.38108. Thus “both forces cooperate for all metrics” would be
incorrect. Both structural recruitment and large-component assembly must be
described separately. Force toggles are model interventions; they do not prove
that one molecular binding pathway exclusively maps to one simulated force.

The `frame_metrics.csv` and four-shell sensitivity table permit checking
assembly over time and dependence on the graph threshold. The simulation ends
at 28 ms. No dense/dilute coexistence, equilibrium stationarity, direct
microscopic contact kinetics, or thermodynamic LLPS has been established by
these measurements alone. Fixed van der Waals attraction is not independently
identified by the factorial design. Report this as a controlled finite-time
organization response of the specified model.
