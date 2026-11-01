# Plotting entry points

Start with [the plotting guide](../docs/plotting.md). Copy and edit
[`configs/plots_9x275.example.json`](../configs/plots_9x275.example.json) once;
all notebooks and command-line renderers read that shared profile.

- `01_reconstruction.ipynb`: mass fit and 1D efficiency.
- `02_w_cross_section.ipynb`: saved W inference.
- `03_differential.ipynb`: saved pT²/rapidity inference and explicit unavailable cases.
- `06_mass_and_2d_efficiency.ipynb`: mass plus two saved efficiency maps.
- `07_mass_and_efficiency_w_xv.ipynb`: latest mass-table uncertainties and W–xV map.
- `render_all.py`: command-line reconstruction/W/differential rendering.

Every render creates fresh outputs. Unvalidated inference requires explicit
`allow_diagnostic_results` or `--allow-diagnostics`. Existing figures and
analysis results are never selected by an automatic latest-run search.

[Source notes](SOURCE_NOTES.md) record which active files were consolidated and
which configuration changes were made. Scientific formulas and figure styles
are preserved.
