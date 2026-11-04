# CP2K Utils

Python tools for analysing and plotting [CP2K](https://www.cp2k.org/) output, focused on point-defect studies.

**What it does**

| | Module | Main entry points |
|---|---|---|
| Charge transition level (CTL) diagrams for crystalline hosts | `ctl` | `plot_ctl` |
| CTL spreads over many sites/cells (amorphous materials) | `ctl` | `plot_range_ctl` |
| Band-aligned CTL summary over several defects and materials | `ctl` | `plot_multi_mat_ctl` |
| Load defect data from a CSV | `ctl` | `load_defects_csv`, `write_sample_csv` |
| Projected density of states (PDOS) and IPR plots | `new_pdos` (current), `pdos` (legacy) | `PDOSPlotter`, `smeared_pdos` |
| XYZ → CIF conversion | `structures` | `xyz_to_cif`, CLI `xyz-to-cif` |
| Charge-density-based interstitial site generation (WIP) | `defect_gen` | `int_gen` |

Worked example notebooks are in [`example_scripts/`](example_scripts/).

---

## Installation

The modules are imported as `cp2k_utils.<module>`, so the folder **that contains** this repository has to be on your Python path:

```bash
git clone https://github.com/isaacmackley/cp2k_utils.git
```

```bash
export PYTHONPATH="$PWD:$PYTHONPATH"   # run from the folder you cloned into
```

Dependencies: `numpy`, `matplotlib`, `pymatgen` (only for `structures` and `defect_gen`). Plots use Arial if it is installed.

---

## Charge transition levels (`cp2k_utils.ctl`)

### Data model

Every input is a small dataclass:

| Class | Represents | Key fields |
|---|---|---|
| `Defect` | one charge state of one defect | `name, charge, energy (Ha), mu_added (Ha), correction (eV), defect_type ('i'/'v'/'s'), site, mu_removed (Ha), label`, optional `title, cell_id, site_id, material` |
| `Cell` | one calculation cell | `cell_id, vbm, cbm, bulk` (Ha), optional `alignment` (Ha) |
| `Material` | a host material made of one or more `Cell`s | `name, title, cells, offset` (conduction band offset vs the reference material, eV) |
| `CTL` | a transition level from outside the package | `transition, energy` (eV above VBM), optional `e_min, e_max` |
| `FermiLevel` | a horizontal Fermi level guide | `label, energy, colour` |
| `StyleOptions` | shared plot styling for every CTL function | fonts, limits, legend, band shading, labels… |

All charge states with the same `label` are one defect. In amorphous data, each (`material`, `label`, `cell_id`, `site_id`) is one configuration. Each configuration is computed against its own cell's bulk energy and VBM.

Formation energies are calculated as

$$E_f(E_F) = E_\mathrm{defect} - E_\mathrm{bulk} - \mu_\mathrm{added} + \mu_\mathrm{removed} + E_\mathrm{corr} + q\,(E_\mathrm{VBM} + E_F)$$

Transition levels are the exact crossings of the lowest-energy lines. When a charge state is skipped (negative-U), the transitions through the skipped state are also reported as `bypassed`.

### `plot_ctl` — single host

```python
from cp2k_utils.ctl import Defect, StyleOptions, plot_ctl

N = -19.82216306984949 / 2       # chemical potential (Ha)
defects = [
    Defect('N',  1, -2614.87186971514484, N, 0.05162, 'i', '', 0, label='N_i'),
    Defect('N',  0, -2614.52029755960666, N, 0.0,     'i', '', 0, label='N_i'),
    Defect('N', -1, -2614.17935614371982, N, 0.05162, 'i', '', 0, label='N_i'),
    # ...
]
results = plot_ctl(defects, bulk=-2604.768017, vbm=0.287532, cbm=0.437028,
                   style=StyleOptions(show_transition_labels=True),
                   save=True, save_format='png', export_csv=False)

results['N_i']['transitions']   # [{'transition': '(1+/1-)', 'E_fermi': 1.60, 'E_form': 4.23, 'kind': 'thermo', ...}, ...]
```

Makes one figure per defect `label`, saved as `'<label> CTL Diagram.<fmt>'`. Each line is solid where that charge state is the most stable and dashed elsewhere. Transition levels are marked, and the VB and CB are shaded. With `export_csv=True` every formation-energy line is also written to a CSV.

[Example notebook →](example_scripts/plot_ctl.ipynb)

<img src="https://github.com/user-attachments/assets/b8155f76-6105-4ee0-b746-360a00581a3a" width="500">
<img src="https://github.com/user-attachments/assets/608a9488-498b-4f81-8239-d9b6616b1f32" width="500">

### `plot_range_ctl` — one defect over many sites

```python
from cp2k_utils.ctl import Cell, Material, load_defects_csv, plot_range_ctl

hfo2 = Material('HfO2', cells=[Cell('a1', 0.2875, 0.4370, -2604.768),
                               Cell('a2', 0.2851, 0.4362, -2604.731)])
defects = load_defects_csv('amorphous_defects.csv')
n_i = [d for d in defects if d.material == 'HfO2' and d.label == 'N_i']

fig, ax = plot_range_ctl(n_i, hfo2, shade=True, title='N$_i$ in a-HfO$_2$')
```

Plots the lowest-energy envelope of every site/cell, shades between the lowest and highest, and draws the overall lowest in bold. The Fermi axis is referenced to the mean VBM of the cells.

[Example notebook →](example_scripts/plot_range_ctl.ipynb)

<img src="https://github.com/user-attachments/assets/5bd64f82-41bc-4b61-ae8d-1c6ffb07c65a" width="500">
<img src="https://github.com/user-attachments/assets/861bcd9c-1ea5-4f9a-8f33-0861db92a75c" width="500">

### `plot_multi_mat_ctl` — band-aligned summary

```python
from cp2k_utils.ctl import Material, FermiLevel, StyleOptions, plot_multi_mat_ctl

materials = [Material('HfO2', cells=hfo2_cells, offset=0.0),     # reference: VBM drawn at 0
             Material('SiO2', cells=sio2_cells, offset=1.5)]     # CBM 1.5 eV above HfO2's

fig, ax, distributions = plot_multi_mat_ctl(
    defects, materials,
    fermi_lines=[FermiLevel('Si E$_F$', 4.2, colour='tab:blue')],
    style=StyleOptions(legend=True, show_bypassed=True),
    save=True, filename='band_alignment.png')
```

- Each material is a block of defect columns, and each block's VB/CB shading is placed by its band offset.
- Each transition is a bar coloured by charge state. With several sites the bar shows min / mean / max.
- Bypassed (negative-U) transitions can be drawn faintly. Labels are spread out automatically so they don't overlap.
- Transition levels calculated elsewhere can be added with `precomputed_ctls={material: {defect: [CTL(...)]}}`.
- `distributions` returns the raw per-site values.

[Example notebook →](example_scripts/plot_multi_mat_ctl.ipynb)

<img src="https://github.com/user-attachments/assets/68b13add-87d2-44b8-9ee9-36a1d3d30c74" width="500">

### Styling

`StyleOptions` is passed as `style=` to every CTL plotting function. Some common options:

```python
StyleOptions(
    figsize=(16, 12), title=None, show_title=False,
    axes_fontsize=30, ticks_fontsize=20, band_alpha=0.6,
    xlim=None, ylim=None,                                    # plot_ctl
    show_transitions=True, show_transition_labels=False,      # plot_ctl
    legend=False, legend_loc='upper left', legend_fontsize=20,
    annotate_ctls=True, ctl_linewidth=3, ctl_colour=None,     # plot_multi_mat_ctl
    show_bypassed=True, bypassed_alpha=0.45,
    show_material_labels=None, band_edge_lines=True,
)
```

See the notebooks for which option applies to which plot.

### CSV input

```python
from cp2k_utils.ctl import write_sample_csv, load_defects_csv
write_sample_csv('template.csv')          # header + example rows
defects = load_defects_csv('my_defects.csv')
```

Required columns: `name, charge, energy, mu_added, correction, defect_type, site, mu_removed, label`.
Optional columns: `title, cell_id, site_id, material`.

---

## PDOS (`cp2k_utils.new_pdos`)

```python
from cp2k_utils.new_pdos import PDOSPlotter

dos = PDOSPlotter(['Hf', 'O', 'N'],      # in CP2K KIND order (k1, k2, ...)
                  spin=True, sigma=0.003, vbm=0, xmin=-3, xmax=6,
                  figw=15, figh=10, output='DOS.png')
dos.plot_individual()   # per element, s/p/d/f resolved   -> <element>_DOS.png
dos.plot_total()        # all elements                   -> DOS.png
dos.plot_ipr()          # total DOS + IPR sticks          -> DOS_IPR.png
```

Run it in the folder with the CP2K `.pdos` files. It matches `*ALPHA_k<n>*` / `*BETA_k<n>*` for spin-polarised runs and `*k<n>*.pdos` otherwise. IPR data are read from `IPR_output-Alpha.dat` / `IPR_output-Beta.dat`. Energies are shifted so the Fermi level is at 0 eV. `smeared_pdos(files, sigma)` returns the smeared data without plotting.

The older function interface `pdos_plot(elements, ...)` / `ipr_pdos_plot(...)` is still available in `cp2k_utils.pdos`.

[Example notebook →](example_scripts/pdos_plot.ipynb)

<img src="https://github.com/user-attachments/assets/b8846064-8d98-492b-b174-e6b52ad99934" width="500">
<img src="https://github.com/user-attachments/assets/7943ff75-a1ac-4440-a306-0e7797927ac6" width="500">

---

## XYZ to CIF (`cp2k_utils.structures`)

Converts an XYZ file to a CIF, given the cell parameters (Å, degrees).

```python
from cp2k_utils.structures import xyz_to_cif
xyz_to_cif('input.xyz', a, b, c, alpha, beta, gamma, 'output.cif')
```

Or from the command line (after installing the console scripts):

```bash
xyz-to-cif input.xyz a b c alpha beta gamma output.cif
```

---

## Interstitial site generator (`cp2k_utils.defect_gen`, WIP)

`int_gen(cube_file, specie, output_filename, int_no=10)` uses pymatgen's `ChargeInterstitialGenerator` on a CP2K density cube file. It writes one CIF per candidate interstitial site, plus one CIF with every site.

---

## Repository layout

```
ctl.py            CTL dataclasses, calculation and plotting
new_pdos.py       PDOSPlotter (current PDOS code)
pdos.py           legacy PDOS functions
structures.py     xyz_to_cif
defect_gen.py     interstitial generator (WIP)
Elements.py       element colour table used in PDOS plots
tools.py          small helpers
cli/              command-line wrappers
example_scripts/  example notebooks and example data
```
