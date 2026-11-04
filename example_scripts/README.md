# Example scripts

Jupyter notebooks showing how to use each part of the package, with an explanation of every input and option.

| Notebook | Covers |
|---|---|
| [plot_ctl.ipynb](./plot_ctl.ipynb) | `Defect` inputs, `plot_ctl`, reading transition levels, `StyleOptions`, CSV input |
| [plot_range_ctl.ipynb](./plot_range_ctl.ipynb) | `Cell` / `Material`, many sites and cells of one defect with `plot_range_ctl` |
| [plot_multi_mat_ctl.ipynb](./plot_multi_mat_ctl.ipynb) | `plot_multi_mat_ctl` with single CTL values (crystalline) and ranges over many sites (amorphous); bar / violin / strip distribution styles, band offsets, Fermi lines, styling, adding precomputed `CTL`s, band-offset-only diagrams |
| [pdos_plot.ipynb](./pdos_plot.ipynb) | `PDOSPlotter` (per-element, total and IPR plots), `smeared_pdos`, legacy `pdos_plot` |

All notebooks are saved with their outputs, so the plots can be viewed on GitHub.

## Example data

**All example data are synthetic.** None of it comes from real calculations. The material names (HfO<sub>2</sub>, SiO<sub>2</sub>) are only there to make the examples readable. Because the numbers are invented, some plots may look odd.

- `data/make_example_data.py` writes `crystalline_defects.csv`, `amorphous_defects.csv` and `amorphous_cells.csv` from made-up formation energies. For the amorphous sets it adds random scatter per cell and per site.
- `data/make_example_pdos.py` writes CP2K-format `.pdos` and IPR files to `data/pdos/`, imitating an N-doped oxide. The plots and smeared data the notebook writes there are git-ignored.
- The band offsets in section 5 of the multi-material notebook and the precomputed SiO<sub>2</sub> levels are also made-up values.

## Running

Make sure the folder that contains `cp2k_utils` is on your `PYTHONPATH` (see the main [README](../README.md)), then run from this directory:

```bash
jupyter notebook
```
