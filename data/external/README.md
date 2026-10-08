# Experimental context and comparison boundary

The relevant experimental article is Hoffmann et al. (2023),
"Synapsin condensation controls synaptic vesicle sequestering and dynamics,"
[Nature Communications](https://doi.org/10.1038/s41467-023-42372-6).
Its author-provided source data are available at
[Figshare DOI 10.6084/m9.figshare.24179445.v1](https://doi.org/10.6084/m9.figshare.24179445.v1)
and from the [publisher's supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-42372-6/MediaObjects/41467_2023_42372_MOESM8_ESM.xlsx).
The unchanged downloaded workbook had SHA-256
`a10e182b2e5eb5781f550cac99c5dde556917b96fe2c198f5b429324cba51f57`.

`hoffmann_mobility_reference.json` contains the extracted mobility
scale-reference values used by the comparison figure, their source and
comparison limitations. The publisher PDF and full third-party workbook are
not redistributed here. The original source provider's terms apply to the
experimental values; the authors' data license does not relicense them.

The study supports the biological relevance of synapsin condensation and
vesicle organization. Its experimental observation windows and geometry do not
match the model's 28-ms trajectories. The comparison code explicitly projects
the model to xy and fits short-lag apparent mobility using slope/4, but this
does not remove the other experimental/model differences. No parameter in the
registered 2 x 2 experiment was retrospectively fitted to the workbook.

The displayed reference is an unmatched scale diagnostic, not an independent
validation or a claim of numerical agreement. Absolute model mobility is
higher than the experimental reference. The simulation outcomes are presented
as finite-time structural model responses, with this mismatch retained and
documented rather than hidden.
