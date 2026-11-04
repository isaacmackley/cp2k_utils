"""
Generates *illustrative* CP2K-format PDOS and IPR files for pdos_plot.ipynb.

These are NOT real calculations. They imitate a spin-polarised, N-doped HfO2
cell: an O-p valence band, a Hf-d conduction band starting ~5.7 eV above the
VBM, and a localised N-p defect state in the gap (occupied in the alpha
channel, empty in the beta channel). Run `python make_example_pdos.py` from
this folder to regenerate the files in ./pdos/.
"""
import os
import numpy as np

rng = np.random.default_rng(1)
Ha = 27.211386
EF_HA = 0.2875                      # Fermi level written in the header (Ha) = VBM
OUT = 'pdos'
os.makedirs(OUT, exist_ok=True)


def ev(x):
    """eV relative to the VBM -> absolute Ha."""
    return EF_HA + np.asarray(x) / Ha


# eigenvalues (eV relative to VBM)
vb = np.sort(-6.5 * rng.beta(1.6, 1.2, 380))          # valence band, -6.5..0
cb = np.sort(5.7 + 6.0 * rng.beta(1.3, 2.0, 300))     # conduction band, 5.7..11.7
gap = {'ALPHA': np.array([0.9]),                      # occupied N-p state
       'BETA': np.array([2.3])}                       # empty N-p state

# (kind index, element, orbitals, weight function -> projections per state)
def projections(element, region, n):
    r = rng.random
    if element == 'Hf':      # s p d f
        if region == 'vb':  w = [0.02, 0.03, 0.25, 0.01]
        elif region == 'cb': w = [0.03, 0.04, 0.70, 0.02]
        else:               w = [0.0, 0.0, 0.05, 0.0]
    elif element == 'O':     # s p
        if region == 'vb':  w = [0.03, 0.65]
        elif region == 'cb': w = [0.03, 0.15]
        else:               w = [0.0, 0.10]
    else:                    # N: s p
        if region == 'vb':  w = [0.0, 0.006]
        elif region == 'cb': w = [0.0, 0.004]
        else:               w = [0.02, 0.80]
    return np.array(w) * (0.5 + r((n, len(w))))


KINDS = [(1, 'Hf', ['s', 'p', 'd', 'f']), (2, 'O', ['s', 'p']), (3, 'N', ['s', 'p'])]

for spin in ('ALPHA', 'BETA'):
    e_all = np.concatenate([vb, gap[spin], cb])
    regions = ['vb'] * len(vb) + ['gap'] * len(gap[spin]) + ['cb'] * len(cb)
    occ = (e_all <= 1.0).astype(float)
    for k, el, orbs in KINDS:
        rows = []
        for reg in ('vb', 'gap', 'cb'):
            n = regions.count(reg)
            rows.append(projections(el, reg, n))
        proj = np.vstack(rows)
        fname = f'{OUT}/HfO2_N-{spin}_k{k}-1_0.pdos'
        with open(fname, 'w') as fh:
            fh.write(f'# Projected DOS for atomic kind {el} at iteration step i = 0, '
                     f'E(Fermi) =     {EF_HA:.6f} a.u.\n')
            fh.write('#     MO Eigenvalue [a.u.]      Occupation ' +
                     ''.join(f'{o:>12}' for o in orbs) + '\n')
            for i, (e, o, p) in enumerate(zip(ev(e_all), occ, proj), 1):
                fh.write(f'{i:8d}{e:16.6f}{o:12.6f}' + ''.join(f'{x:12.5f}' for x in p) + '\n')

    # IPR: small for delocalised band states, large for the defect state
    ipr = np.where(np.array(regions) == 'gap', 0.35, 0.005 + 0.02 * rng.random(len(e_all)))
    name = 'Alpha' if spin == 'ALPHA' else 'Beta'
    np.savetxt(f'{OUT}/IPR_output-{name}.dat', np.column_stack([ev(e_all), ipr]), fmt='%.6f')

print('wrote', sorted(os.listdir(OUT)))
