# Example scripts

Jupyter notebooks showing how to use each part of the package, with an explanation of every input and option.

| Notebook | Covers |
|---|---|
| [plot_ctl.ipynb](./plot_ctl.ipynb) | `Defect` inputs, `plot_ctl`, reading transition levels, `StyleOptions`, CSV input |
| [plot_range_ctl.ipynb](./plot_range_ctl.ipynb) | `Cell` / `Material`, many sites and cells of one defect with `plot_range_ctl` |
| [plot_multi_mat_ctl.ipynb](./plot_multi_mat_ctl.ipynb) | Band offsets, `plot_multi_mat_ctl`, Fermi lines, styling, adding precomputed `CTL`s |
| [pdos_plot.ipynb](./pdos_plot.ipynb) | `PDOSPlotter` (per-element, total and IPR plots), `smeared_pdos`, legacy `pdos_plot` |

The CTL notebooks are saved with their outputs, so the plots can be viewed on GitHub. The PDOS notebook needs CP2K `.pdos` files and is saved without outputs.

## Example data

- `data/amorphous_defects.csv` and `data/amorphous_cells.csv` are used by the range and multi-material notebooks. They are **illustrative only**: `data/make_example_data.py` made them by adding random noise to the crystalline HfO<sub>2</sub> values in `plot_ctl.ipynb`. They are not real calculations.

## Running

Make sure the folder that contains `cp2k_utils` is on your `PYTHONPATH` (see the main [README](../README.md)), then run from this directory:

```bash
jupyter notebook
```
