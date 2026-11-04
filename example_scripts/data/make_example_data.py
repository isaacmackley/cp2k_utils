"""
Generates the SYNTHETIC defect data used by the CTL example notebooks.

Nothing here comes from a real calculation. Each defect is defined by made-up
target formation energies (eV, at E_F = VBM) for each charge state; these are
converted into the 'total energies' that the package expects, and for the
amorphous sets random noise is added per cell and per site. Material names
(HfO2, SiO2) are only there to make the examples readable.

Run `python make_example_data.py` in this folder to regenerate:
    crystalline_defects.csv   one cell, one site per defect
    amorphous_defects.csv     several cells x sites per defect, two materials
    amorphous_cells.csv       band edges / bulk energy of each amorphous cell
"""
import csv
import numpy as np

Ha = 27.211386
rng = np.random.default_rng(7)

MU = {'N': -10.0, 'O': -16.0}          # chemical potentials (Ha)


def correction(q):
    return round(0.06 * q ** 2, 5)     # made-up charge correction (eV)


# defect key -> (name, type, site, mu_added, mu_removed, label, title)
DEFECTS = {
    'N_i': ('N', 'i', '',  MU['N'], 0.0,     'N_i', ''),
    'N_O': ('N', 's', 'O', MU['N'], MU['O'], 'N_O', ''),
    'V_O': ('O', 'v', 'O', 0.0,     MU['O'], 'V_O', ''),
}

# target formation energy at the VBM (eV) per charge state
TARGETS = {
    'HfO2': {
        'N_i': {2: 1.0, 1: 1.8, 0: 3.4, -1: 4.6, -2: 8.5},    # negative-U (1+/1-)
        'N_O': {2: 3.5, 1: 3.8, 0: 4.5, -1: 6.0, -2: 9.5},
        'V_O': {2: -1.0, 1: 2.5, 0: 5.6, -1: 10.4, -2: 15.6},  # negative-U (2+/0)
    },
    'SiO2': {
        'N_i': {2: 0.5, 1: 2.0, 0: 5.0, -1: 8.5},
        'V_O': {2: 1.0, 1: 4.0, 0: 7.5, -1: 13.0},
    },
}

HEADER = ['name', 'charge', 'energy', 'mu_added', 'correction', 'defect_type',
          'site', 'mu_removed', 'label', 'title', 'cell_id', 'site_id', 'material']


def total_energy(ef_vbm, q, bulk, vbm, mu_added, mu_removed):
    """Invert E_f(VBM) = (E - bulk - mu_a + mu_r)*Ha + corr + q*VBM*Ha."""
    return bulk + mu_added - mu_removed + (ef_vbm - correction(q) - q * vbm * Ha) / Ha


def row(key, q, energy, cell_id, site_id, material):
    name, dtype, site, mu_a, mu_r, label, title = DEFECTS[key]
    return [name, q, round(energy, 8), mu_a, correction(q), dtype, site, mu_r,
            label, title, cell_id, site_id, material]


# ---------------- crystalline: one cell, exact targets ----------------
CRYSTAL = dict(vbm=0.25, cbm=0.46, bulk=-1000.0)       # Ha (gap ~5.7 eV)
rows = []
for key, targets in TARGETS['HfO2'].items():
    for q, ef in targets.items():
        e = total_energy(ef, q, CRYSTAL['bulk'], CRYSTAL['vbm'], *DEFECTS[key][3:5])
        rows.append(row(key, q, e, '', '', 'c-HfO2'))
with open('crystalline_defects.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(HEADER)
    w.writerows(rows)

# ---------------- amorphous: several cells x sites, with noise ----------------
AMORPHOUS = {   # material -> (n_cells, mean vbm, mean cbm, bulk, sites per cell)
    'HfO2': (3, 0.25, 0.455, -1000.0, 8),
    'SiO2': (2, 0.20, 0.530, -700.0, 8),
}
EDGE_SCATTER = 0.25 / Ha          # cell-to-cell band edge scatter (Ha)
SITE_SCATTER = 0.25               # site-to-site shift of all charge states (eV)
CHARGE_SCATTER = 0.12             # extra scatter per charge state (eV)

defect_rows, cell_rows = [], []
for mat, (n_cells, vbm0, cbm0, bulk0, n_sites) in AMORPHOUS.items():
    for c in range(n_cells):
        cell_id = f'{mat}_c{c + 1}'
        vbm = round(vbm0 + rng.normal(0, EDGE_SCATTER), 6)
        cbm = round(cbm0 + rng.normal(0, EDGE_SCATTER), 6)
        bulk = round(bulk0 + rng.normal(0, 0.05), 6)
        cell_rows.append([mat, cell_id, vbm, cbm, bulk])
        for key, targets in TARGETS[mat].items():
            for s in range(n_sites):
                site_id = f'{cell_id}_s{s + 1}'
                shift = rng.normal(0, SITE_SCATTER)
                for q, ef in targets.items():
                    ef_site = ef + shift + rng.normal(0, CHARGE_SCATTER)
                    e = total_energy(ef_site, q, bulk, vbm, *DEFECTS[key][3:5])
                    defect_rows.append(row(key, q, e, cell_id, site_id, mat))

with open('amorphous_defects.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(HEADER)
    w.writerows(defect_rows)
with open('amorphous_cells.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['material', 'cell_id', 'vbm', 'cbm', 'bulk'])
    w.writerows(cell_rows)

print(f'crystalline: {len(rows)} rows; amorphous: {len(defect_rows)} rows, {len(cell_rows)} cells')
