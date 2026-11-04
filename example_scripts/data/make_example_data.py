"""
Generates the *illustrative* amorphous-ensemble CSV used by the
plot_range_ctl and plot_multi_mat_ctl example notebooks.

The numbers are NOT real calculations. They are the crystalline HfO2 N_i and
V_O energies from plot_ctl.ipynb, copied onto several made-up cells/sites with
random noise added, so the example plots have a realistic-looking spread.
Re-run with `python make_example_data.py` to regenerate the files.
"""
import csv
import numpy as np

rng = np.random.default_rng(42)

N = -19.82216306984949 / 2      # chemical potential of N (Ha)
O = -31.80451126833290 / 2      # chemical potential of O (Ha)

REF_BULK = -2604.76801711694679
REF_VBM = 0.287532

# (name, defect_type, site, mu_added, mu_removed, label, title)
DEFECTS = {
    'N_i': ('N', 'i', '', N, 0.0, 'N_i', None),
    'V_O': ('O', 'v', 'O', 0.0, O, 'V_O', None),
}
# reference energies for each charge state (Ha) and charge corrections (eV)
REF_ENERGIES = {
    'N_i': {2: (-2615.18759465001131, 0.20648), 1: (-2614.87186971514484, 0.05162),
            0: (-2614.52029755960666, 0.0), -1: (-2614.17935614371982, 0.05162)},
    'V_O': {2: (-2589.39963503920580, 0.20648), 1: (-2589.00866444766143, 0.05162),
            0: (-2588.61414089194932, 0.0), -1: (-2588.18012071979138, 0.05162)},
}

# material -> list of (cell_id, vbm, cbm, bulk) in Ha
MATERIALS = {
    'HfO2': [('a1', 0.2875, 0.4370, -2604.768),
             ('a2', 0.2851, 0.4362, -2604.731),
             ('a3', 0.2899, 0.4381, -2604.802)],
    'SiO2': [('s1', 0.2010, 0.5330, -1530.412),
             ('s2', 0.1985, 0.5302, -1530.388)],
}
SITES_PER_CELL = 3

defect_rows, cell_rows = [], []
for mat, cells in MATERIALS.items():
    for cell_id, vbm, cbm, bulk in cells:
        cell_rows.append([mat, cell_id, vbm, cbm, bulk])
        for key, (name, dtype, site, mu_a, mu_r, label, title) in DEFECTS.items():
            for s in range(SITES_PER_CELL):
                site_id = f'{cell_id}_site{s + 1}'
                site_shift = rng.normal(0, 0.01)          # site-to-site scatter (Ha)
                for q, (e_ref, corr) in REF_ENERGIES[key].items():
                    # keep formation energies/CTLs similar to the reference,
                    # re-referenced to this cell's bulk energy and VBM
                    e = (e_ref - REF_BULK + bulk - q * (vbm - REF_VBM)
                         + site_shift + rng.normal(0, 0.004))
                    defect_rows.append([name, q, round(e, 8), mu_a, corr, dtype,
                                        site, mu_r, label, title or '',
                                        cell_id, site_id, mat])

with open('amorphous_defects.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['name', 'charge', 'energy', 'mu_added', 'correction',
                'defect_type', 'site', 'mu_removed', 'label', 'title',
                'cell_id', 'site_id', 'material'])
    w.writerows(defect_rows)

with open('amorphous_cells.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['material', 'cell_id', 'vbm', 'cbm', 'bulk'])
    w.writerows(cell_rows)

print(f'wrote {len(defect_rows)} defect rows and {len(cell_rows)} cells')
