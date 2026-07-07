import re
import csv
import warnings

from dataclasses import dataclass, field
from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator
from collections import defaultdict


plt.rcParams['font.family'] = 'Arial'
Ha_to_eV = 27.211386


##### COLOURS, TRANSITIONS, HELPER FUNCTIONS #####

# colors = plt.cm.tab10.colors
# colors = ['#0072B2','#E69F00','#009E73','#D55E00','#CC79A7','#56B4E9','#F0E442']
colors = ['tab:blue', 'tab:red', 'tab:cyan', 'tab:pink', 'tab:olive', 'tab:brown', 'tab:purple', 'tab:gray']

CHARGE_COLOURS = {
    4: '#00B3B3',   # teal
    3: '#0F2080',   # navy
    2: '#1E88E5',   # blue
    1: '#7A1F1F',   # dark red / maroon
    0: '#000000',   # black (neutral - reads as the reference state)
   -1: '#C71585',   # magenta
   -2: '#785EF0',   # violet
   -3: '#8C6D1F',   # dark gold / ochre
   -4: '#00B3B3',   # teal
}

FALLBACK_COLOURS = ['#0F2080', '#1E88E5', '#7A1F1F', '#C71585',
                    '#785EF0', '#8C6D1F', '#00B3B3', '#000000']

def charge_colour(charge: int) -> str:
    """Stable colour for a given charge state."""
    if charge in CHARGE_COLOURS:
        return CHARGE_COLOURS[charge]
    return FALLBACK_COLOURS[abs(charge) % len(FALLBACK_COLOURS)]

def charge_label(charge: int) -> str:
    """LaTeX charge label: '$1^+$', '$0$', '$2^-$'."""
    if charge > 0:
        return f'${abs(charge)}^+$'
    if charge < 0:
        return f'${abs(charge)}^-$'
    return '$0$'

def _parse_charge(tok: str) -> int:
    """'2+' -> 2, '+1' -> 1, '0' -> 0, '3-' -> -3, '-2' -> -2. Tolerates
    surrounding parentheses/whitespace from compute_ctls labels like '(+1/0)'."""
    tok = tok.strip().strip('()').strip()
    if tok in ('0', '+0', '-0'):
        return 0
    m = re.fullmatch(r'([+-]?)(\d+)([+-]?)', tok)
    if not m:
        raise ValueError(f'Cannot parse charge token {tok!r}')
    pre, num, post = m.groups()
    sign = pre or post or '+'
    return int(num) * (-1 if sign == '-' else 1)

def _parse_transition(s: str) -> tuple:
    """'2+/1+' -> (2, 1); '0/-1' -> (0, -1)."""
    a, b = s.split('/')
    return _parse_charge(a), _parse_charge(b)
 
 
def transition_colour(transition: str) -> str:
    """Colour for a transition, keyed on the higher charge of the pair."""
    return charge_colour(max(_parse_transition(transition)))

def _format_species(name):
    """Render a species name with its digits as subscripts: 'N2' -> 'N$_2$'."""
    return re.sub(r'(\d+)', r'$_{\1}$', name)


def _has_digit(name):
    return any(ch.isdigit() for ch in name)


def format_defect_label(name, defect_type, site='', title=None):
    """Build a display label (LaTeX) from the raw fields. `title` overrides all.
        'N','i'           -> 'N$_{i}$'
        'N2','i'          -> '(N$_{2}$)$_{i}$'
        'C','s',site='N'  -> 'C$_{N}$'
        'O','v'           -> 'V$_{O}$'
    """
    if title:
        return title
    t = (defect_type or '').lower()
    species = _format_species(name)
    wrapped = f'({species})' if _has_digit(name) else species
    if t == 'i':
        return f'{wrapped}$_{{i}}$'
    if t == 's':
        sub = site if site else 's'
        return f'{wrapped}$_{{{sub}}}$'
    if t == 'v':
        return f'V$_{{{wrapped}}}$'
    return species


##### DATACLASSES #####

_DEFAULT_CELL = 'cell0'
_DEFAULT_SITE = 'site0'

@dataclass
class Defect:
    name: str
    charge: int
    energy: float       # Ha
    mu_added: float     # Ha
    correction: float   # eV
    defect_type: str    # 'i', 's', 'v'
    site: str
    mu_removed: float   # Ha
    label: str
    title: Optional[str] = None
    # stuff for amorphous
    cell_id: str = _DEFAULT_CELL
    site_id: str = _DEFAULT_SITE
    material: Optional[str] = None

    @property
    def key(self) -> str:
        """Single source of truth for the per-charge dictionary key."""
        return f'{self.name}_{self.defect_type}{self.charge}'
    
    def base_formation_energy(self, bulk: float, e_vbm: float) -> float:
        """
        Charge-independent part of the formation energy, in eV. Full curve:
            E_form(E_fermi) = base + charge * (e_vbm + E_fermi)
        All inputs in Ha except `correction` (eV) and `e_vbm` (eV).
        """
        return (
            self.energy * Ha_to_eV
            - bulk * Ha_to_eV
            - self.mu_added * Ha_to_eV
            + self.mu_removed * Ha_to_eV
            + self.correction
        )


@dataclass
class Cell:
    """info for one calculation cell"""
    cell_id: str
    vbm: float              # Ha
    cbm: float              # Ha
    bulk: float             # Ha
    alignment: float = 0.0  # Ha

    @property
    def vbm_eV(self) -> float:
        return (self.vbm + self.alignment) * Ha_to_eV
 
    @property
    def cbm_eV(self) -> float:
        return (self.cbm + self.alignment) * Ha_to_eV
 
    @property
    def bulk_eV(self) -> float:
        return self.bulk * Ha_to_eV



@dataclass
class Material:
    name: str
    title: Optional[str] = None
    cells: list = field(default_factory=list)
    offset: float = 0.0     # eV


    def cell(self, cell_id: str) -> Cell:
        if len(self.cells) == 1:
            return self.cells[0]
        for c in self.cells:
            if c.cell_id == cell_id:
                return c
        raise KeyError(f"Material '{self.name}' has no cell '{cell_id}'")
 
    @property
    def display_label(self) -> str:
        """How the material is shown on the plot. `name` stays the plain ASCII
        identifier (dict key, reference= argument); this is the rendered form."""
        return self.title if self.title else _format_species(self.name)

    @property
    def mean_vbm_eV(self) -> float:
        return float(np.mean([c.vbm_eV for c in self.cells]))
 
    @property
    def mean_cbm_eV(self) -> float:
        return float(np.mean([c.cbm_eV for c in self.cells]))
 
    @property
    def mean_gap_eV(self) -> float:
        return self.mean_cbm_eV - self.mean_vbm_eV
 
    @property
    def vbm_spread_eV(self) -> tuple:
        vs = [c.vbm_eV for c in self.cells]
        return min(vs), max(vs)
 
    @property
    def cbm_spread_eV(self) -> tuple:
        cs = [c.cbm_eV for c in self.cells]
        return min(cs), max(cs)
    

    def cbm_on_axis(self, ref_cbm_disp: float) -> float:
        return ref_cbm_disp + self.offset
 
    def vbm_on_axis(self, ref_cbm_disp: float) -> float:
        return self.cbm_on_axis(ref_cbm_disp) - self.mean_gap_eV
 
    def _axis_shift(self, ref_cbm_disp: float) -> float:
        """Maps this material's raw eV scale onto the displayed axis so its mean
        CBM lands at cbm_on_axis (preserving per-cell spread)."""
        return self.cbm_on_axis(ref_cbm_disp) - self.mean_cbm_eV
 
    def vbm_spread_on_axis(self, ref_cbm_disp: float) -> tuple:
        s = self._axis_shift(ref_cbm_disp)
        lo, hi = self.vbm_spread_eV
        return lo + s, hi + s
 
    def cbm_spread_on_axis(self, ref_cbm_disp: float) -> tuple:
        s = self._axis_shift(ref_cbm_disp)
        lo, hi = self.cbm_spread_eV
        return lo + s, hi + s
 
    def ctl_on_axis(self, ctl_value_rel_vbm: float, ref_cbm_disp: float) -> float:
        """CTL (eV above this material's VBM) -> displayed-axis position."""
        return self.vbm_on_axis(ref_cbm_disp) + ctl_value_rel_vbm



@dataclass
class CTL:
    """A computed (or externally supplied) charge transition level, in eV
    relative to the host VBM."""
    transition: str
    energy: float
    e_min: Optional[float] = None
    e_max: Optional[float] = None
 
    @property
    def has_spread(self) -> bool:
        return self.e_min is not None and self.e_max is not None
 
    @property
    def colour(self) -> str:
        return transition_colour(self.transition)


@dataclass
class FermiLevel:
    label: str
    energy: float
    colour: Optional[str] = None


@dataclass
class StyleOptions:

    # --- figure / axes ---
    figsize: tuple = None            # None -> per-function default
    title: str = None                # overrides the function's default title
    show_title: bool = False         # plot_ctl: draw the title at all
    band_alpha: float = 0.6          # VB/CB shading opacity
    axes_fontsize: float = 30
    ticks_fontsize: float = 20
 
    # --- CTL marks / lines ---
    ctl_linewidth: float = 3
    ctl_label_fontsize: float = 20
    annotate_ctls: bool = True      # label each CTL next to its mark
 
    # --- legend ---
    legend: bool = False
    legend_loc: str = 'upper left'   # plot_ctl
    legend_frameon: bool = False     # plot_ctl
    legend_fontsize: float = 20
 
    # --- plot_ctl transitions / limits ---
    show_transitions: bool = True
    show_transition_labels: bool = False
    xlim: tuple = None
    ylim: tuple = None
 
    # --- multi-mat specifics ---
    band_edge_lines: bool = True
    show_material_labels: bool = None   # None -> auto (only if >1 material)
    mats_fontsize: float = 30
    show_bypassed: bool = True
    bypassed_linewidth: float = None    # None -> match ctl_linewidth
    bypassed_alpha: float = 0.45
    leader_linewidth: float = 1.0



##### Other useful functions etc. #####



def gradient_fill(ax, x_start, x_end, ymin, ymax, color, fade_direction, alpha_max=0.3):
    """
    Adds a horizontal gradient fill between x_start and x_end.
    fade_direction: 'left'  — opaque at x_start, fades right
                    'right' — opaque at x_end,   fades left
    """
    rgba = np.zeros((1, 256, 4))
    base = mcolors.to_rgba(color)
    
    alphas = np.linspace(alpha_max, 0, 256) if fade_direction == 'left' else np.linspace(0, alpha_max, 256)
    
    rgba[0, :, :3] = base[:3]  # RGB constant across gradient
    rgba[0, :, 3]  = alphas    # alpha varies

    ax.imshow(rgba, aspect='auto', extent=[x_start, x_end, ymin, ymax],
              origin='lower', zorder=0, interpolation='bilinear')

def _vgradient(ax, x_start, x_end, y_lo, y_hi, color, fade_direction,
               alpha_max=0.3):
    """Vertical gradient fill between y_lo and y_hi. fade_direction: 'up' ->
    opaque at y_lo fading up; 'down' -> opaque at y_hi fading down."""
    if y_hi <= y_lo:
        return
    rgba = np.zeros((256, 1, 4))
    rgba[:, 0, :3] = mcolors.to_rgba(color)[:3]
    col = (np.linspace(alpha_max, 0, 256) if fade_direction == 'up'
           else np.linspace(0, alpha_max, 256))
    rgba[:, 0, 3] = col
    ax.imshow(rgba, aspect='auto', extent=[x_start, x_end, y_lo, y_hi],
              origin='lower', zorder=0, interpolation='bilinear')

def _group_configurations(defects):
    """
    Bucket a flat Defect list into configurations.
    Returns dict keyed by (material, label, cell_id, site_id) -> list[Defect].
    A configuration is one site in one cell across its charge states.
    """
    groups = defaultdict(list)
    for d in defects:
        groups[(d.material, d.label, d.cell_id, d.site_id)].append(d)
    return groups
 
def _pair_crossing(by_charge, qa, qb, e_vbm):
    """Analytic Fermi-level crossing of charge states qa and qb (relative to
    VBM), and the formation energy there. Returns (ef_cross, e_form)."""
    ba, bb = by_charge[qa], by_charge[qb]
    ef_cross = (bb - ba) / (qa - qb) - e_vbm
    e_form = ba + qa * (e_vbm + ef_cross)
    return float(ef_cross), float(e_form)

def compute_ctls(config_defects, e_vbm, bulk):
    """
    Compute charge transition levels for ONE configuration (one site, one cell).
 
    Parameters
    ----------
    config_defects : list[Defect]   all charge states of one site/cell
    e_vbm : float                   that cell's VBM in eV (incl. alignment)
    bulk : float                    that cell's pristine energy in Ha
 
    Returns
    -------
    list[dict]  {'transition': '(+1/0)', 'E_fermi': float, 'E_form': float}
                E_fermi is relative to the cell VBM.
    """
    by_charge = {}
    for d in config_defects:
        base = d.base_formation_energy(bulk, e_vbm)
        by_charge[d.charge] = base   # slope is charge; line: base + q*(e_vbm+Ef)
 
    charges = sorted(by_charge.keys(), reverse=True)  # high -> low charge
    if len(charges) < 2:
        return []
 
    Ef = np.linspace(-20.0, 20.0, 40001)
 
    def line(q):
        return by_charge[q] + q * (e_vbm + Ef)
 
    stacked = np.vstack([line(q) for q in charges])
    lowest = np.array([charges[i] for i in np.argmin(stacked, axis=0)])
    change = np.where(np.diff(lowest) != 0)[0]
 
    transitions = []
    for idx in change:
        q1 = int(lowest[idx])
        q2 = int(lowest[idx + 1])
        ef_cross, e_form = _pair_crossing(by_charge, q1, q2, e_vbm)

        lo_q, hi_q = sorted((q1, q2))
        bypassed = [q for q in charges if lo_q < q < hi_q]

        transitions.append({
            'transition': f'({q1:+}/{q2:+})',
            'E_fermi': ef_cross,
            'E_form': e_form,
            'kind': 'thermo',
            'skipped': bool(bypassed),
            'bypassed': bypassed,
        })

        if bypassed:
            seq = [q1] + sorted(bypassed, reverse=(q1 > q2)) + [q2]
            for qa, qb in zip(seq[:-1], seq[1:]):
                ef_b, ef_form_b = _pair_crossing(by_charge, qa, qb, e_vbm)
                transitions.append({
                    'transition': f'({qa:+}/{qb:+})',
                    'E_fermi': ef_b,
                    'E_form': ef_form_b,
                    'kind': 'bypassed',
                    'skipped': False,
                    'bypassed': [],
                })

    return transitions

def summarize(values):
    """min/mean/max of a list of values -> (lo, mean, hi)."""
    arr = np.asarray(values, dtype=float)
    return float(arr.min()), float(arr.mean()), float(arr.max())
 
 
def set_integer_yticks(ax, ymin, ymax, fontsize=20):
    """Robust integer y-ticks."""
    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(axis='y', labelsize=fontsize)



##### CSV stuff #####



_CSV_FIELDS = {
    'name': str, 'charge': int, 'energy': float, 'mu_added': float,
    'correction': float, 'defect_type': str, 'site': str,
    'mu_removed': float, 'label': str, 'title': str,
    'cell_id': str, 'site_id': str, 'material': str,
}
_CSV_REQUIRED = ['name', 'charge', 'energy', 'mu_added', 'correction',
                 'defect_type', 'site', 'mu_removed', 'label']
 
 
def load_defects_csv(path):
    """
    Load a flat list[Defect] from a CSV. See _CSV_FIELDS for the layout; a
    header row is required. Optional columns fall back to dataclass defaults.
    """
    defects = []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in _CSV_REQUIRED if c not in reader.fieldnames]
        if missing:
            raise ValueError(f'CSV missing required columns: {missing}')
        for row in reader:
            kwargs = {}
            for col, caster in _CSV_FIELDS.items():
                if col not in row or row[col] is None:
                    continue
                val = row[col]
                if val == '':
                    # Keep empty strings for str fields (e.g. an interstitial's
                    # blank `site`); skip empties for numeric/optional fields so
                    # they fall back to dataclass defaults.
                    if caster is str and col in _CSV_REQUIRED:
                        kwargs[col] = ''
                    continue
                kwargs[col] = caster(val)
            defects.append(Defect(**kwargs))
    return defects
 
 
def write_sample_csv(path):
    """Write a small example CSV illustrating the expected layout."""
    header = list(_CSV_FIELDS.keys())
    rows = [
        # name,charge,energy,mu_added,correction,defect_type,site,mu_removed,
        # label,title,cell_id,site_id,material
        ['N', 1, -100.50, 0.5, 0.10, 'i', '', 0.0, 'N_i', 'N$_i$',
         'cellA', 'siteA', 'HfO2'],
        ['N', 0, -100.80, 0.5, 0.00, 'i', '', 0.0, 'N_i', 'N$_i$',
         'cellA', 'siteA', 'HfO2'],
        ['N', -1, -101.00, 0.5, 0.10, 'i', '', 0.0, 'N_i', 'N$_i$',
         'cellA', 'siteA', 'HfO2'],
    ]
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)





##### Main Plotting Functions #####


def plot_ctl(defects, bulk, vbm, cbm, style=None, save=True, show=False, save_format='pdf', export_csv=False):
    """
    Plots CTL diagrams for input defects.

    defects: list of Defect dataclass instances
    bulk: energy of pristine cell (Ha)
    vbm: VBM energy (Ha)
    cbm: CBM energy (Ha)
    """
    
    opts = style or StyleOptions()
    figsize = opts.figsize or (15, 10)

    E_cbm = cbm * Ha_to_eV
    E_vbm = vbm * Ha_to_eV
    band_gap = E_cbm - E_vbm
    PAD= 1
    
    xmax = round(band_gap + PAD)
    E_fermi = np.linspace(-PAD, xmax, ((xmax + PAD) * 100) + 1)


    groups = defaultdict(list)
    for d in defects:
        groups[d.label].append(d)
    
    results = {}


    for group_label, group in groups.items():
        first = group[0]
 
        E_curves = {}
        E_low = None
        ymin_c, ymax_c = [], []
        for d in group:
            base = d.base_formation_energy(bulk, E_vbm)
            curve = base + d.charge * (E_vbm + E_fermi)
            E_curves[d.key] = curve
            E_low = curve if E_low is None else np.minimum(E_low, curve)
            ymin_c += [base + d.charge * E_vbm, base + d.charge * E_cbm]
            ymax_c.append(base + d.charge * E_vbm)
 
        ymin = round(min(ymin_c) - 1)
        ymax = round(max(ymax_c) + 1)
 
        title_text = (first.title or first.label) + ' Charge Levels'
 
        fig, ax = plt.subplots(figsize=(15, 10))
        if opts.show_title:
            ax.set_title(title_text, size=40, pad=30)
        ax.set_xlabel('Fermi Energy (eV)', size=opts.axes_fontsize, labelpad=30)
        ax.set_ylabel('Formation Energy (eV)', size=opts.axes_fontsize, labelpad=30)
 
        ax.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
        ax.axvline(x=band_gap, color='tab:orange', linestyle='-', alpha=0.5)
        gradient_fill(ax, -PAD, 0, ymin, ymax, 'tab:green',
                      fade_direction='right', alpha_max=opts.band_alpha)
        gradient_fill(ax, band_gap, xmax, ymin, ymax, 'tab:orange',
                      fade_direction='left', alpha_max=opts.band_alpha)
 
        for d in group:
            curve = E_curves[d.key]
            colour = charge_colour(d.charge)
            is_lowest = np.isclose(curve, E_low, atol=1e-6)
            ax.plot(E_fermi, np.where(is_lowest, curve, np.nan),
                    lw=3, color=colour, label=charge_label(d.charge), zorder=4)
            ax.plot(E_fermi, np.where(~is_lowest, curve, np.nan),
                    lw=3, color=colour, linestyle='--', alpha=0.7, zorder=3)
 
        # Transitions via the shared per-configuration routine.
        transitions = compute_ctls(group, E_vbm, bulk)
        if opts.show_transitions:
            for t in transitions:
                if t['kind'] != 'thermo':
                    continue
                if -PAD <= t['E_fermi'] <= band_gap + PAD:
                    ax.axvline(x=t['E_fermi'], color='k', linestyle=':',
                               lw=1, alpha=0.4, zorder=2)
                    ax.plot(t['E_fermi'], t['E_form'], 'ko',
                            markersize=10, zorder=5)
                    if opts.show_transition_labels:
                        ax.annotate(rf"$\epsilon${t['transition']}",
                                    xy=(t['E_fermi'], ymin + 0.2),
                                    ha='center', fontsize=opts.ctl_label_fontsize)
 
        ax.set_xlim(opts.xlim if opts.xlim is not None else (-PAD, xmax))
        ax.tick_params(axis='x', labelsize=opts.ticks_fontsize)
        if opts.ylim is not None:
            ax.set_ylim(opts.ylim)
            ax.tick_params(axis='y', labelsize=opts.ticks_fontsize)
        else:
            set_integer_yticks(ax, ymin, ymax, fontsize=opts.ticks_fontsize)
 
        ax.tick_params(axis='both', which='major', direction='in',
                       length=6, width=1.5, labelsize=opts.ticks_fontsize)
        ax.tick_params(axis='both', which='minor', direction='in', length=3)
        ax.minorticks_on()
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        ax.legend(frameon=opts.legend_frameon, markerscale=5.0,
                  fontsize=opts.legend_fontsize, loc=opts.legend_loc)
        ax.grid(visible=False)
 
        if save:
            fig.savefig(f'{first.label} CTL Diagram.{save_format}',
                        bbox_inches='tight', dpi=300)
        if export_csv:
            header = ['E_fermi'] + [charge_label(d.charge).strip('$')
                                    for d in group]
            data = np.column_stack([E_fermi] + [E_curves[d.key] for d in group])
            np.savetxt(f'{first.label} CTL data.csv', data, delimiter=',',
                       header=','.join(header), comments='')
        if show:
            plt.show()
 
        results[group_label] = {'figure': fig, 'transitions': transitions}
 
    return results


    # type_labels = {
    #     'i': lambda d: f"{d.name}_i",
    #     's': lambda d: f"{d.name}_{d.site}",
    #     'v': lambda d: f"V_{d.site}",
    # }

    # types = defaultdict(list)
    # for defect in defects:
    #     types[defect.label].append(defect)

    # E_low = {}
    # defect_types = []

    # for type_key, typegroup in types.items():

    #     first = typegroup[0]
    #     defect_type = type_labels[first.defect_type](first)

    #     E_ = {}
    #     ymin1 = float('inf')
    #     ymin2 = float('inf')
    #     ymax_group = float('-inf')

    #     for i, defect in enumerate(typegroup):
    #         name = defect.name + '_' + defect.defect_type + str(defect.charge)
    #         base_E = (defect.energy * Ha_to_eV - bulk * Ha_to_eV - defect.mu_added * Ha_to_eV + defect.mu_removed * Ha_to_eV + defect.correction)

    #         E_[name] = base_E + defect.charge * (E_vbm + E_fermi)

    #         if i == 0:
    #             defect_types.append(defect_type)
    #             E_low[defect_type] = E_[name]
    #         else:
    #             E_low[defect_type] = np.minimum(E_low[defect_type], E_[name])

    #         e_at_vbm = base_E + defect.charge * E_vbm
    #         e_at_cbm = base_E + defect.charge * E_cbm

    #         ymin1 = min(ymin1, e_at_vbm)
    #         ymax_group = max(ymax_group, e_at_vbm)
    #         ymin2 = min(ymin2, e_at_cbm)

    #     ymin = round(min(ymin1, ymin2) - 1)
    #     ymax = round(ymax_group + 1)

    #     title = str(type_labels[typegroup[0].defect_type](typegroup[0])) + ' Charge Levels'

    #     generated_label = type_labels[first.defect_type](first)
    #     plot_key = first.label if first.label is not None else generated_label
    #     opts = plot_options.get(plot_key, default_options)

    #     plt.figure(figsize=(15, 10))
    #     if title == True:
    #         plt.title(title, size=40, pad=30)
    #     plt.xlabel('Fermi Energy (eV)', size=30, labelpad=30)
    #     plt.ylabel('Formation Energy (eV)', size=30, labelpad=30)

    #     ax = plt.gca()

    #     plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
    #     plt.axvline(x=E_cbm - E_vbm, color='tab:orange', linestyle='-', alpha=0.5)
    #     gradient_fill(ax, -1, 0, ymin, ymax, 'tab:green', fade_direction='right', alpha_max=0.6)
    #     gradient_fill(ax, E_cbm - E_vbm, xmax, ymin, ymax, 'tab:orange', fade_direction='left', alpha_max=0.6)
        
    #     for i, defect in enumerate(typegroup):
    #         name = defect.name + '_' + defect.defect_type + str(defect.charge)
    #         color = colors[i % len(colors)]
            
    #         if defect.charge > 0:
    #             label = f'${abs(defect.charge)}^+$'
    #         elif defect.charge < 0:
    #             label = f'${abs(defect.charge)}^-$'
    #         else:
    #             label = 'neutral'
            
    #         is_lowest = np.isclose(E_[name], E_low[defect_type], atol=1e-6)
    #         solid  = np.where(is_lowest,  E_[name], np.nan)
    #         dashed = np.where(~is_lowest, E_[name], np.nan)

    #         plt.plot(E_fermi, solid,  marker='', lw=3, color=color, label=label, zorder=4)
    #         plt.plot(E_fermi, dashed, marker='', lw=3, color=color, linestyle='--', alpha=0.7, zorder=3)

    #     transitions = []

    #     # Find which charge state is lowest at each Fermi energy point
    #     lowest_charge_at_each_point = []
    #     for ef_idx in range(len(E_fermi)):
    #         lowest_charge = min(
    #             typegroup,
    #             key=lambda d: E_[d.name + '_' + d.defect_type + str(d.charge)][ef_idx]
    #         ).charge
    #         lowest_charge_at_each_point.append(lowest_charge)

    #     lowest_charge_at_each_point = np.array(lowest_charge_at_each_point)

    #     # Find where the lowest charge state changes
    #     change_indices = np.where(np.diff(lowest_charge_at_each_point) != 0)[0]

    #     for idx in change_indices:
    #         q1 = lowest_charge_at_each_point[idx]
    #         q2 = lowest_charge_at_each_point[idx + 1]

    #         name1 = typegroup[0].name + '_' + typegroup[0].defect_type + str(q1)
    #         name2 = typegroup[0].name + '_' + typegroup[0].defect_type + str(q2)

    #         # Precise crossing via linear interpolation
    #         diff = E_[name1] - E_[name2]
    #         x1, x2 = E_fermi[idx], E_fermi[idx + 1]
    #         y1, y2 = diff[idx], diff[idx + 1]
    #         E_cross = x1 - y1 * (x2 - x1) / (y2 - y1)
    #         E_form_cross = np.interp(E_cross, E_fermi, E_[name1])

    #         transitions.append({
    #             'E_fermi': E_cross,
    #             'E_form': E_form_cross,
    #             'label': f'$\\epsilon$({q1:+}/{q2:+})'
    #         })

        

    #     band_gap = E_cbm - E_vbm

    #     for trans in transitions:
    #         # Only show transitions within the band gap
    #         if -1 <= trans['E_fermi'] <= band_gap + 1:
    #             ax.axvline(x=trans['E_fermi'], color='k', linestyle=':', lw=1, alpha=0.4, zorder=2)
    #             ax.plot(trans['E_fermi'], trans['E_form'], 'ko', markersize=10, zorder=5)
    #             # ax.annotate(
    #             # trans['label'],
    #             # xy=(trans['E_fermi'], ymin + 0.2),
    #             # ha='center',
    #             # fontsize=14
    #             # )

    #     if opts.xlim is not None:
    #         plt.xlim(opts.xlim)
    #     else:
    #         plt.xlim([-1, xmax])

    #     plt.xticks(fontsize=20)

    #     if opts.ylim is not None:
    #         plt.ylim(opts.ylim)
    #         plt.yticks(fontsize=20)
    #     else:
    #         plt.ylim([ymin, ymax])
    #         if ymin == 0 or ymax == 0:
    #             plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) + 1), fontsize=20)
    #         else:
    #             if ymin / abs(ymin) == ymax / abs(ymax):
    #                 plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) - 1), fontsize=20)
    #             else:
    #                 plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) + 1), fontsize=20)

    #     plt.tick_params(axis='both', which='major', direction='in', length=6, width=1.5, labelsize=20)
    #     plt.tick_params(axis='both', which='minor', direction='in', length=3)
    #     plt.minorticks_on()

        
    #     for spine in ax.spines.values():
    #         spine.set_linewidth(1.5)

    #     plt.legend(frameon=opts.legend_frameon, markerscale=5.0, fontsize=20, loc=opts.legend_loc)
    #     plt.grid(visible=0)

    #     plt.savefig(f'{title} CTL Diagram.pdf', bbox_inches='tight', dpi=300)


def plot_range_ctl(defects, material, shade=True, title=None,
                   save=True, show=False, save_format='pdf', filename=None):
    """
    Overlay formation-energy lines for every site/cell of one defect type
    (one amorphous ensemble). The spread across sites is the visual; CTL dots
    are intentionally omitted (a crossing is only defined per-site).
 
    Parameters
    ----------
    defects : list[Defect]
        All charge states of all sites/cells of ONE defect type. Each Defect's
        cell_id must resolve to a Cell in `material`.
    material : Material
        Provides per-cell vbm/cbm/bulk. Curves are referenced to the mean VBM
        across the cells so sites sit on a common axis (per-site physics stays
        exact; only the shared zero is averaged).
    shade : bool
        Fill the band between the lowest and highest formation energy across
        sites at each Fermi level.
    """
    configs = _group_configurations(defects)
    if not configs:
        raise ValueError('No defects supplied to plot_range_ctl.')
 
    mean_vbm = material.mean_vbm_eV
    mean_cbm = material.mean_cbm_eV
    band_gap = mean_cbm - mean_vbm
    PAD = 1.0
    xmax = round(band_gap + PAD)
    E_fermi = np.linspace(-PAD, xmax, int((xmax + PAD) * 100) + 1)
 
    # Each configuration -> its ground-state envelope, referenced to mean VBM.
    envelopes = []
    label = None
    for (mat, lab, cell_id, site_id), rows in configs.items():
        label = label or lab
        cell = material.cell(cell_id)
        e_vbm = cell.vbm_eV
        env = None
        for d in rows:
            base = d.base_formation_energy(cell.bulk, e_vbm)
            # Reference Fermi axis to the mean VBM: replace e_vbm with mean_vbm
            # in the (e_vbm + Ef) term so all configs share one zero, while the
            # base (site physics) stays exact.
            curve = base + d.charge * (mean_vbm + E_fermi)
            env = curve if env is None else np.minimum(env, curve)
        envelopes.append(env)
 
    envelopes = np.vstack(envelopes)
    lo = envelopes.min(axis=0)
    hi = envelopes.max(axis=0)
 
    ymin = round(float(envelopes.min()) - 1)
    ymax = round(float(envelopes.max()) + 1)
 
    fig, ax = plt.subplots(figsize=(15, 10))
    if title:
        ax.set_title(title, size=40, pad=30)
    ax.set_xlabel('Fermi Energy (eV)', size=30, labelpad=30)
    ax.set_ylabel('Formation Energy (eV)', size=30, labelpad=30)
 
    # VB/CB guides at the mean edges.
    ax.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
    ax.axvline(x=band_gap, color='tab:orange', linestyle='-', alpha=0.5)
    gradient_fill(ax, -PAD, 0, ymin, ymax, 'tab:green',
                  fade_direction='right', alpha_max=0.3)
    gradient_fill(ax, band_gap, xmax, ymin, ymax, 'tab:orange',
                  fade_direction='left', alpha_max=0.3)
 
    colour = '#1E88E5'
    if shade:
        ax.fill_between(E_fermi, lo, hi, color=colour, alpha=0.2, lw=0, zorder=2)
    for env in envelopes:
        ax.plot(E_fermi, env, lw=1.5, color=colour, alpha=0.6, zorder=3)
    # Bold the overall lowest envelope (the physically observed ground state).
    ax.plot(E_fermi, lo, lw=3, color=colour, zorder=4)
 
    ax.set_xlim(-PAD, xmax)
    set_integer_yticks(ax, ymin, ymax, fontsize=20)
    ax.tick_params(axis='both', which='major', direction='in',
                   length=6, width=1.5, labelsize=20)
    ax.tick_params(axis='both', which='minor', direction='in', length=3)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.grid(visible=False)
 
    if filename is None:
        filename = f'{label} Range CTL Diagram.{save_format}'
    if save:
        fig.savefig(filename, bbox_inches='tight', dpi=300)
    if show:
        plt.show()
 
    return fig, ax


# def plot_range_ctl(defects,bulks,legend=False):
#     """
#     Plots CTL range diagrams for input defects

#     defects: list of lists with format ['Defect Specie', Charge, Defect Cell Energy, Chemical Potential of added species, Err correction, Defect type ('i', 's', 'v' for interstial, substional or vacancy), Site name (e.g. O for an O site. Can be left as '' for interstial), Chemical potential for removed species)
#     bulk: Energy of pristine cell
#     vmb: vbm energy
#     cbm: cbm energy

#     note - all energies should be given in Ha
#     """

#     Ha_to_eV = 27.211386

#     E_vbm = 0
#     E_cbm = 0
#     n=0
#     while n < len(bulks):
#         E_vbm = E_vbm + (bulks[n][2] * Ha_to_eV)
#         E_cbm = E_cbm + (bulks[n][3] * Ha_to_eV)
#         n=n+1
#     E_vbm_ave = E_vbm/len(bulks)
#     E_cbm_ave = E_cbm/len(bulks)
#     xmax=round((E_cbm_ave-E_vbm_ave)+1)
#     E_fermi = np.linspace(-1,xmax,((xmax+1)*100)+1)


#     types = {}
#     n=0
#     while n < len(defects):
#         type_key = defects[n][2]+defects[n][7]+defects[n][8]

#         if type_key not in types:
#             types[type_key] = []
#         types[type_key].append(defects[n])

#         n=n+1
    
#     E_low={}
#     E_lo={}
#     E_hi={}
#     defect_types = []
#     defect_charges = []

#     for type_key, typegroup in types.items():

        
    
#         E_={}
#         xequal0_={}
#         xequalcbm_={}
#         ymin1 = float(1000)
#         ymin2 = float(1000)
#         ymax = float(-1000)

#         m=0

        

#         for defect in typegroup:   
            
#             n=0
#             while n<len(bulks):
#                 if bulks[n][0] == defect[0]:
#                     bulk=bulks[n][1]
#                     vbm=bulks[n][2]*Ha_to_eV
#                     cbm=bulks[n][3]*Ha_to_eV
#                     n=len(bulks)
#                 n=n+1
            
#             name = 'cell'+str(defect[0])+'site'+str(defect[1])+defect[2]+'_'+defect[7]+str(defect[3])
#             E_[name] = defect[4]*Ha_to_eV - bulk*Ha_to_eV - defect[5]*Ha_to_eV + defect[9]*Ha_to_eV + defect[3] * (vbm + E_fermi) + defect[6]
            
#             if typegroup[0][7] == 'i':
#                 defect_type = str(defect[2])+'_'+str(defect[7])
#             if typegroup[0][7] == 's':
#                 defect_type = str(defect[2])+'_'+str(defect[8])
#             if typegroup[0][7] == 'v':
#                 defect_type = 'V_'+str(defect[8])
            
            
#             defect_charge = defect_type+str(defect[3])
            
#             if defect_charge not in defect_charges:
#                 defect_charges.append(defect_charge)
#             if defect_charge not in E_lo:
#                 E_lo[defect_charge] = E_[name]
#             else:
#                 E_lo[defect_charge] = np.minimum(E_lo[defect_charge], E_[name])
#             if defect_charge not in E_hi:
#                 E_hi[defect_charge] = E_[name]
#             else:
#                 E_hi[defect_charge] = np.maximum(E_hi[defect_charge], E_[name])

#             if m == 0:
#                 defect_types.append(defect_type)
#                 E_low[defect_type] = E_[name]
#             else:
#                 E_low[defect_type] = np.minimum(E_low[defect_type], E_[name])


#             xequal0_[name] = defect[4]*Ha_to_eV - bulk*Ha_to_eV - defect[5]*Ha_to_eV + defect[9]*Ha_to_eV + defect[3] * (vbm) + defect[6]
#             ymin1 = np.minimum(ymin1, xequal0_[name])
#             ymax = np.maximum(ymax, xequal0_[name])

#             xequalcbm_[name] = defect[4]*Ha_to_eV - bulk*Ha_to_eV - defect[5]*Ha_to_eV + defect[9]*Ha_to_eV + defect[3] * (cbm) + defect[6]
#             ymin2 = np.minimum(ymin2, xequalcbm_[name])

#             m=m+1
           

#         ymin = np.minimum(ymin1, ymin2)
#         ymin = round(ymin - 1)
#         ymax = round(ymax + 1)

#         plt.figure(figsize=(15,10))
#         if typegroup[0][7] == 'i':
#             plt.title(f'{typegroup[0][2]}$_{'i'}$ Charge Levels', size=40, pad=30)
#         if typegroup[0][7] == 's':
#             plt.title(f'{typegroup[0][2]}$_{{{typegroup[0][8]}}}$ Charge Levels', size=40, pad=30)
#         if typegroup[0][7] == 'v':
#             plt.title(f'V$_{{{typegroup[0][8]}}}$ Charge Levels', size=40, pad=30)
#         plt.xlabel('Fermi Energy (eV)', size=30, labelpad=30) 
#         plt.ylabel('Formation Energy (eV)', size=30, labelpad=30)

#         #VBM and CBM 
#         plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
#         plt.axvline(x=E_cbm_ave-E_vbm_ave, color='tab:orange', linestyle='-', alpha=0.5)
#         plt.fill([-1,0,0,-1], [ymin,ymin,ymax,ymax], 'tab:green',[E_cbm_ave-E_vbm_ave,xmax,xmax,E_cbm_ave-E_vbm_ave], [ymin,ymin,ymax,ymax], 'tab:orange', alpha=0.2)

#         #Defect Plots
#         for defect in typegroup:
            
#             name = 'cell'+str(defect[0])+'site'+str(defect[1])+defect[2]+'_'+defect[7]+str(defect[3])
            
#             q=-5
#             while q < 6:
#                 color = colors[q % len(colors)]
#                 if defect[3] == q:
#                     if defect[3] == 0:
#                         plt.plot(E_fermi, E_[name], marker = '', lw=3, color=color, label='neutral')
#                     elif defect[3] > 0:
#                         plt.plot(E_fermi, E_[name], marker = '', lw=3, color=color, label=f'{str(abs(defect[3]))}{'+'}')
#                     elif defect[3] < 0:
#                         plt.plot(E_fermi, E_[name], marker = '', lw=3, color=color, label=f'{str(abs(defect[3]))}{'-'}')
#                 q=q+1


#         #Lowest state line
#         #plt.plot(E_fermi, x_int, linestyle='--', color='k', lw=2, label='')

#         #Other details
#         plt.xlim([-1,xmax])
#         plt.ylim([ymin,ymax])
#         plt.xticks(np.linspace(-1,xmax,abs(xmax)+2),fontsize=20)
#         if ymin == 0 or ymax == 0:
#             plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1),fontsize=20)
#         else:
#             if ymin/abs(ymin) == ymax/abs(ymax):
#                 plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)-1),fontsize=20)
#             else:
#                 plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1),fontsize=20)

#         plt.grid(visible=0)#, which='major', axis='both',linestyle='--')

#         #if legend==True:
#             #plt.legend(markerscale=5.0, fontsize=15, loc='upper left')

#         if typegroup[0][7] == 'i':
#             plt.savefig(f'{typegroup[0][2]}_i Multi CTL Diagram', bbox_inches='tight')
#         if typegroup[0][7] == 's':
#             plt.savefig(f'{typegroup[0][2]}_{typegroup[0][8]} Multi CTL Diagram', bbox_inches='tight')
#         if typegroup[0][7] == 'v':
#             plt.savefig(f'V_{typegroup[0][8]} Multi CTL Diagram', bbox_inches='tight')
    



#         plt.figure(figsize=(15,10))
#         if typegroup[0][7] == 'i':
#             plt.title(f'{typegroup[0][2]}$_{'i'}$ Charge Levels', size=30, pad=30)
#         if typegroup[0][7] == 's':
#             plt.title(f'{typegroup[0][2]}$_{{{typegroup[0][8]}}}$ Charge Levels', size=30, pad=30)
#         if typegroup[0][7] == 'v':
#             plt.title(f'V$_{{{typegroup[0][8]}}}$ Charge Levels', size=30, pad=30)
#         plt.xlabel('Fermi Energy (eV)', size=20, labelpad=30) 
#         plt.ylabel('Formation Energy (eV)', size=20, labelpad=30)

#         #VBM and CBM 
#         plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
#         plt.axvline(x=E_cbm_ave-E_vbm_ave, color='tab:orange', linestyle='-', alpha=0.5)
#         plt.fill([-1,0,0,-1], [ymin,ymin,ymax,ymax], 'tab:green',[E_cbm_ave-E_vbm_ave,xmax,xmax,E_cbm_ave-E_vbm_ave], [ymin,ymin,ymax,ymax], 'tab:orange', alpha=0.2)

#         #Defect Plots
#         for defect in typegroup:

#             name = 'cell'+str(defect[0])+'site'+str(defect[1])+defect[2]+'_'+defect[7]+str(defect[3])
            
#             n=0
#             while n < len(defect_charges):
#                 color = colors[n % len(colors)]
#                 plt.plot(E_fermi, E_lo[defect_charges[n]], marker = '', lw=1, color=color, label=f'{defect_charges[n]}')
#                 plt.plot(E_fermi, E_hi[defect_charges[n]], marker = '', lw=1, color=color)

#                 plt.fill([-1,xmax,xmax,-1], [np.interp(-1, E_fermi, E_lo[defect_charges[n]]),np.interp(xmax, E_fermi, E_lo[defect_charges[n]]),np.interp(xmax, E_fermi, E_hi[defect_charges[n]]),np.interp(-1, E_fermi, E_hi[defect_charges[n]])], color=color, alpha=0.2)

#                 n=n+1

#         #Lowest state line
#         #plt.plot(E_fermi, x_int, linestyle='--', color='k', lw=2, label='')

#         #Other details
#         plt.xlim([-1,xmax])
#         plt.ylim([ymin,ymax])
#         plt.xticks(np.linspace(-1,xmax,abs(xmax)+2))
#         if ymin == 0 or ymax == 0:
#             plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1))
#         else:
#             if ymin/abs(ymin) == ymax/abs(ymax):
#                 plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)-1))
#             else:
#                 plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1))

#         plt.grid(visible=0)#, which='major', axis='both',linestyle='--')

#         #if legend==True:
#             #plt.legend(markerscale=5.0, fontsize=15, loc='upper left')

#         if typegroup[0][7] == 'i':
#             plt.savefig(f'{typegroup[0][2]}_i Range CTL Diagram', bbox_inches='tight')
#         if typegroup[0][7] == 's':
#             plt.savefig(f'{typegroup[0][2]}_{typegroup[0][8]} Range CTL Diagram', bbox_inches='tight')
#         if typegroup[0][7] == 'v':
#             plt.savefig(f'V_{typegroup[0][8]} Range CTL Diagram', bbox_inches='tight')
    

def _ensemble_ctls(defects, material):
    """
    For one defect type in one material, compute per-site CTLs and collect them
    by transition. Returns dict:
        transition -> {'values': [E_fermi...], 'kind': 'thermo'|'bypassed',
                       'skipped': bool}
    All E_fermi relative to each site's own cell VBM.
    """
    configs = _group_configurations(defects)
    agg = {}
    for (mat, lab, cell_id, site_id), rows in configs.items():
        cell = material.cell(cell_id)
        for t in compute_ctls(rows, cell.vbm_eV, cell.bulk):
            tr = t['transition']
            if tr not in agg:
                agg[tr] = {'values': [], 'kind': t['kind'],
                           'skipped': t.get('skipped', False)}
            agg[tr]['values'].append(t['E_fermi'])
            if t['kind'] == 'thermo':
                agg[tr]['kind'] = 'thermo'
                agg[tr]['skipped'] = agg[tr]['skipped'] or t.get('skipped', False)
    return agg


def _place_labels(ax, specs, ymin, ymax, fontsize, x_pad=0.03, leader_lw=1.0):
    """
    Draw CTL annotations for one column, nudging them apart vertically so they
    don't overlap, and connecting each to its line with a thin leader.
    specs: list of {'text','y','colour','alpha','x'}.
    """
    if not specs:
        return
    span = ymax - ymin
    min_gap = 1.6 * (fontsize / 72.0) / max(ax.figure.get_size_inches()[1], 1) * span

    order = sorted(range(len(specs)), key=lambda i: specs[i]['y'])
    ys = [specs[i]['y'] for i in order]

    placed = ys[:]
    for k in range(1, len(placed)):
        if placed[k] - placed[k - 1] < min_gap:
            placed[k] = placed[k - 1] + min_gap
    overshoot = placed[-1] - (ymax - 0.1)
    if overshoot > 0:
        placed = [p - overshoot for p in placed]
    for k in range(len(placed) - 2, -1, -1):
        if placed[k + 1] - placed[k] < min_gap:
            placed[k] = placed[k + 1] - min_gap
    placed = [min(max(p, ymin + 0.1), ymax - 0.1) for p in placed]

    label_x = max(s['x'] for s in specs) + x_pad
    for slot, i in enumerate(order):
        s = specs[i]
        y_line = s['y']
        y_lab = placed[slot]
        if abs(y_lab - y_line) > 1e-3:
            ax.plot([s['x'], label_x], [y_line, y_lab], lw=leader_lw,
                    color=s['colour'], alpha=0.5 * s['alpha'], zorder=3)
        ax.annotate(s['text'], xy=(label_x, y_lab), va='center', ha='left',
                    fontsize=fontsize, color=s['colour'], alpha=s['alpha'])

 
def plot_multi_mat_ctl(defects, materials, fermi_lines=None,
                       style_kind='bar', precomputed_ctls=None, reference=None,
                       style=None, save=False, show=True, filename=None):
    """
    Cross-material CTL summary with band alignment.
 
    For each defect type in each material, per-site CTLs are computed exactly
    (relative to each site's own cell VBM - never mixing charge states across
    sites) and shown as a distribution per transition.
 
    Parameters
    ----------
    defects : list[Defect]
        Flat list across all materials/defects/sites/cells. Each Defect needs
        `material` set (and cell_id/site_id for amorphous ensembles).
    materials : list[Material]
        Each owns its Cells. Band shading uses the min-max spread of cell edges.
    fermi_lines : FermiLevel | list[FermiLevel] | float | None
        Horizontal dashed Fermi guide(s).
    style : 'bar'
        Distribution display. 'bar' draws min/mean/max marks (today). The raw
        per-site values are retained internally, so 'violin'/'strip' can be
        added here later without touching the computation.
    precomputed_ctls : dict | None
        Optional bypass for external CTLs you did not compute yourself:
            {material_name: {defect_name: [CTL, ...]}}
        Used instead of computing from raw Defects for those entries.
    reference : str | None
        Material name whose VBM defines the displayed y = 0 (energy increases
        upward). Defaults to the material with offset == 0 (the CBO reference),
        else the first material. Fermi lines are interpreted on this same axis.
 
    Returns
    -------
    (fig, ax, distributions)
        distributions: {material: {defect: {transition: [values...]}}}
    """

    opts = style or StyleOptions()
    figsize = opts.figsize or (20,15)
    title = opts.title

    if style_kind != 'bar':
        warnings.warn(f"style_kind={style_kind!r} not yet implemented; using 'bar'.",
                      stacklevel=2)
        style_kind = 'bar'
 
    # Normalise fermi_lines.
    if fermi_lines is None:
        fermi_lines = []
    elif isinstance(fermi_lines, FermiLevel):
        fermi_lines = [fermi_lines]
    elif isinstance(fermi_lines, (int, float)):
        fermi_lines = [FermiLevel('E$_F$', float(fermi_lines))]
 
    mat_by_name = {m.name: m for m in materials}
 
    # Displayed axis: reference material's VBM = 0, energy increases upward.
    # The reference is the named material if given, else the one with offset==0,
    # else the first. With ref VBM at 0, the reference CBM sits at its own gap;
    # every other material is placed by its CBO relative to that.
    if reference is not None:
        ref_mat = mat_by_name[reference]
    else:
        zero_off = [m for m in materials if m.offset == 0.0]
        ref_mat = zero_off[0] if zero_off else materials[0]
    ref_cbm_disp = ref_mat.mean_gap_eV   # reference CBM on the displayed axis
 
    # Defect names present per material, in input order.
    # Defect identity for columns is the `label` (e.g. 'V_O'), consistent with
    # plot_ctl grouping - not the bare element `name` ('V').
    defect_order = defaultdict(list)
    label_title = {}
    for d in defects:
        if d.material not in mat_by_name:
            warnings.warn(f"Defect '{d.label}' references unknown material "
                          f"'{d.material}'.", stacklevel=2)
            continue
        if d.label not in defect_order[d.material]:
            defect_order[d.material].append(d.label)
        label_title.setdefault(d.label, format_defect_label(d.name, d.defect_type, d.site, d.title))
    if precomputed_ctls:
        for mat, dd in precomputed_ctls.items():
            for dname in dd:
                if dname not in defect_order[mat]:
                    defect_order[mat].append(dname)
                label_title.setdefault(d.label, format_defect_label(d.name, d.defect_type, d.site, d.title))
 
    # Compute the distribution of CTLs per (material, defect).
    # distributions[material][defect][transition] = [E_fermi values]
    distributions = {}
    for m in materials:
        distributions[m.name] = {}
        for dname in defect_order[m.name]:
            if precomputed_ctls and dname in precomputed_ctls.get(m.name, {}):
                bt = {}
                for c in precomputed_ctls[m.name][dname]:
                    vals = ([c.e_min, c.energy, c.e_max] if c.has_spread
                            else [c.energy])
                    if c.transition not in bt:
                        bt[c.transition] = {'values': [], 'kind': 'thermo',
                                            'skipped': False}
                    bt[c.transition]['values'].extend(vals)
                distributions[m.name][dname] = bt
            else:
                subset = [d for d in defects
                          if d.material == m.name and d.label == dname]
                distributions[m.name][dname] = _ensemble_ctls(subset, m)
 
    # y-range from cell edges across all materials, on the displayed axis.
    all_v = [m.vbm_spread_on_axis(ref_cbm_disp)[0] for m in materials]
    all_c = [m.cbm_spread_on_axis(ref_cbm_disp)[1] for m in materials]
    ctl_axis_vals = []
    for mname, dmap in distributions.items():
        mm = mat_by_name[mname]
        for entry in dmap.values():
            for tr in entry.values():
                for v in tr['values']:
                    ctl_axis_vals.append(mm.ctl_on_axis(v, ref_cbm_disp))
    lo_all = min(all_v + ctl_axis_vals) if ctl_axis_vals else min(all_v)
    hi_all = max(all_c + ctl_axis_vals) if ctl_axis_vals else max(all_c)
    ymin = round(lo_all - 1)
    ymax = round(hi_all + 1)
 
    # x-layout: contiguous block of defect columns per material + separator.
    columns = []
    material_spans = []
    x = 1
    for mi, m in enumerate(materials):
        names = defect_order[m.name]
        start = x - 0.5
        if names:
            for dname in names:
                columns.append((x, dname, m.name))
                x += 1
        else:
            x += 1
        end = x - 0.5
        material_spans.append((m, start, end, (start + end) / 2))
        # if mi != len(materials) - 1:
        #     x += 1            # separator only between materials
    xmin = 0.5
    xmax = x - 0.5
 
    fig, ax = plt.subplots(figsize=figsize)
    if title:
        ax.set_title(title, size=50, pad=40)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel('Fermi Energy (eV)', size=opts.axes_fontsize, labelpad=20)
    ax.set_xlabel('Defects', size=opts.axes_fontsize, labelpad=30)
    ax.set_xticks([c[0] for c in columns])
    ax.set_xticklabels([label_title.get(c[1], c[1]) for c in columns],
                       fontsize=opts.ticks_fontsize)
    ax.tick_params(axis='y', labelsize=opts.ticks_fontsize)

    # Material labels: by default shown only when there is more than one
    # material (a single label is redundant). Override with True/False.
    if opts.show_material_labels is None:
        show_labels = len(materials) > 1
    else:
        show_labels = opts.show_material_labels
 
    # Band shading: min-max spread of cell edges per material, on the common
    # (CBO-aligned) axis. The faded strips show genuine cell-to-cell edge
    # variation; the solid regions are the conservative band extents.
    n_span = len(material_spans)
    for si, (m, start, end, centre) in enumerate(material_spans):
        x0 = xmin if si == 0 else start
        x1 = xmax if si == n_span - 1 else end
        v_lo, v_hi = m.vbm_spread_on_axis(ref_cbm_disp)
        c_lo, c_hi = m.cbm_spread_on_axis(ref_cbm_disp)
        _vgradient(ax, x0, x1, ymin, v_lo, 'tab:green', 'down', opts.band_alpha)
        ax.fill_between([x0, x1], v_lo, v_hi, color='tab:green',
                        alpha=opts.band_alpha * 0.5, lw=0, zorder=0)
        _vgradient(ax, x0, x1, c_hi, ymax, 'tab:orange', 'up', opts.band_alpha)
        ax.fill_between([x0, x1], c_lo, c_hi, color='tab:orange',
                        alpha=opts.band_alpha * 0.5, lw=0, zorder=0)
        if opts.band_edge_lines:
            v_mean = 0.5 * (v_lo + v_hi)
            c_mean = 0.5 * (c_lo + c_hi)
            ax.hlines(v_mean, x0, x1, color='tab:green', linestyle='-',
                      lw=1.5, alpha=0.6, zorder=1)
            ax.hlines(c_mean, x0, x1, color='tab:orange', linestyle='-',
                      lw=1.5, alpha=0.6, zorder=1)
        if show_labels:
            ax.annotate(m.display_label, [centre, ymin + 0.5], fontsize=opts.mats_fontsize,
                        ha='center')
 
    # CTL distribution marks. Each value is a Fermi energy above that material's
    # own VBM; map it onto the common (CBO-aligned) axis before plotting.
    seen_transitions = {}     # legend: thermo transitions only
    for x_centre, dname, mname in columns:
        m = mat_by_name[mname]
        dist = distributions[mname][dname]
        thermo_items = [(tr, d) for tr, d in dist.items()
                        if d.get('kind', 'thermo') == 'thermo']
        bypassed_items = [(tr, d) for tr, d in dist.items()
                          if d.get('kind') == 'bypassed']
        thermo_items.sort(key=lambda kv: max(_parse_transition(kv[0])),
                          reverse=True)
        n_t = len(thermo_items)

        label_specs = []

        for j, (trans, entry) in enumerate(thermo_items):
            values = entry['values']
            if not values:
                continue
            w = 0.0 if n_t <= 1 else 0.4 / (n_t - 1)
            xm = x_centre + w * (j - (n_t - 1) / 2)
            colour = transition_colour(trans)
            seen_transitions.setdefault(trans, colour)

            common_vals = [m.ctl_on_axis(v, ref_cbm_disp) for v in values]
            lo, mean, hi = summarize(common_vals)
            lw = opts.ctl_linewidth
            if lo == hi:
                ax.hlines(mean, xm - 0.1, xm + 0.1, lw=lw, colors=colour,
                          zorder=5)
                right = xm + 0.1
            else:
                ax.fill_between([xm - 0.05, xm + 0.05], lo, hi,
                                color=colour, alpha=0.2, lw=0)
                ax.hlines(lo, xm - 0.075, xm + 0.075, lw=lw, colors=colour, zorder=5)
                ax.hlines(hi, xm - 0.075, xm + 0.075, lw=lw, colors=colour, zorder=5)
                ax.hlines(mean, xm - 0.05, xm + 0.05, lw=lw, colors=colour, zorder=5)
                right = xm + 0.075
            if opts.annotate_ctls:
                label_specs.append({'text': trans, 'y': mean, 'colour': colour,
                                    'alpha': 1.0, 'x': right})

        if opts.show_bypassed:
            for trans, entry in bypassed_items:
                values = entry['values']
                if not values:
                    continue
                colour = transition_colour(trans)
                common_vals = [m.ctl_on_axis(v, ref_cbm_disp) for v in values]
                _, mean, _ = summarize(common_vals)
                blw = opts.bypassed_linewidth if opts.bypassed_linewidth is not None \
                    else opts.ctl_linewidth
                ax.hlines(mean, x_centre - 0.1, x_centre + 0.1,
                          lw=blw, colors=colour,
                          linestyle=(0, (4, 2)), alpha=opts.bypassed_alpha, zorder=4)
                if opts.annotate_ctls:
                    label_specs.append({'text': trans, 'y': mean,
                                        'colour': colour, 'alpha': 0.6,
                                        'x': x_centre + 0.1})

        if opts.annotate_ctls and label_specs:
            _place_labels(ax, label_specs, ymin, ymax,
                          fontsize=opts.ctl_label_fontsize,
                          leader_lw=opts.leader_linewidth)
 
    for fl in fermi_lines:
        ax.axhline(y=fl.energy, label=fl.label, linestyle='--', lw=3,
                   alpha=0.8, color=fl.colour)
 
    if opts.legend:
        handles, labels = [], []
        for trans, colour in sorted(
                seen_transitions.items(),
                key=lambda kv: max(_parse_transition(kv[0])), reverse=True):
            handles.append(Line2D([0], [0], color=colour, lw=6))
            labels.append(trans)
        for fl in fermi_lines:
            handles.append(Line2D([0], [0], color=fl.colour or 'k',
                                   lw=3, linestyle='--'))
            labels.append(fl.label)
        if handles:
            ax.legend(handles, labels, loc='center right', fontsize=opts.legend_fontsize,
                      frameon=True)
 
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(axis='both', which='major', direction='in',
                   length=6, width=1.5, labelsize=opts.ticks_fontsize)
    ax.tick_params(axis='y', which='minor', direction='in', length=3)
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', length=0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.grid(visible=False)
 
    if filename is None:
        filename = 'Multi_Mat_CTLs.png'
    if save:
        fig.savefig(filename, bbox_inches='tight', dpi=300)
    if show:
        plt.show()
 
    return fig, ax, distributions


# def plot_multi_mat_ctl(defects,materials,fermi_line=None):
#     """
#     Plots all given CTLs in given materials with band alignment.

#     nb
#         Ensure that the VBM and CBM energies are given relative to one another between different materials.
#         Ensure that the 'Material Name' is the same in the defects and materials lists.
#         CTL level energies to be given relative to the VBM of the material.

#     args:

#     defects:
#             Each row: [Defect Name, Material Name, List of CTL tuples [Name, Energy] ]
#             e.g.
#                 species = [
#                     [f'N$_{'i'}$', f'HfO$_2$', [
#                         ['+1/0', 1.3469378771661913],
#                         ['0/-1', 2.8837208880316783],
#                         ['-1/-2', 2.369073993248366],
#                         ['-2/-3', 2.6696762876888496],
#                     ]],
#                     [f'N$_{'O'}$', f'HfO$_2$', [
#                         ['3+/2+', 1.7950036300732979],
#                         ['2+/1+', 3.0614261580351125],
#                         ['+1/0', 2.0040361516655647],
#                         ['0/-1', 2.164946988290245],
#                     ]],
#                     [f'N$_{'2'}$$_{'i'}$', f'HfO$_2$', [
#                         ['0/-1', 3.380293837285865],
#                         ['-1/-2', 4.198164515355801],
#                         ['-2/-3', 4.441756959205842],
#                         ['-3/-4', 3.7187798290658316]
#                     ]]
#                 ]
#     materials:
#                 Each row: [Material Name, VBM, CBM]
#                 e.g.
#                     materials = [
#                         ['Si', Sivbm, Sicbm],
#                         [f'HfO$_{'2'}$', HfO2vbm, HfO2cbm]
#                     ]
#     fermi_line:
#                 Draws dashed line at desired fermi level
#     """

#     defects_list = []
#     materials_list = []

#     for row in defects:
#         defects_list.append(row[:2])
#     for i in range(len(materials)):
#         materials_list.append(materials[i][0])

#     ymin = round(min(materials, key=lambda x: x[1])[1]-1)
#     ymax = round(max(materials, key=lambda x: x[2])[2]+1)

#     xmax = len(defects)+len(materials)

#     plt.figure(figsize=(20,20))

#     plt.ylim(ymin,ymax)
#     plt.xlim(0,xmax)

#     n=0
#     for i in range(len(materials)):
#         no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
#         defects_list.insert(n+no_of_defects,['','',''])
#         n=n+no_of_defects+1
#     defects_list.insert(0,['','',''])

#     plt.xticks(np.linspace(0,xmax,xmax+1), [row[0] for row in defects_list], fontsize=20)
#     plt.yticks(fontsize=20)

#     plt.ylabel('Fermi Energy (eV)', size=30, labelpad=30) 

#     n=0
#     for i in range(len(materials)):

#         no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
#         xmin_i = n
#         xmax_i = n+no_of_defects+1
#         vbm_i = materials[i][1]
#         cbm_i = materials[i][2]

#         # plt.hlines(
#         #     y=vbm_i,
#         #     xmin=xmin_i,
#         #     xmax=xmax_i,
#         #     color='tab:green',
#         #     linestyle='-',
#         #     alpha=0.5
#         # )
        
#         # plt.hlines(
#         #     y=cbm_i,
#         #     xmin=xmin_i,
#         #     xmax=xmax_i,
#         #     color='tab:orange',
#         #     linestyle='-',
#         #     alpha=0.5
#         # )

#         plt.fill(
#             [xmin_i,xmin_i,xmax_i,xmax_i],
#             [ymin,vbm_i,vbm_i,ymin],
#             'tab:green',
#             [xmin_i,xmin_i,xmax_i,xmax_i],
#             [cbm_i,ymax,ymax,cbm_i],
#             'tab:orange',
#             alpha=0.2
#         )

#         x_i = n + (no_of_defects+1)/2

#         plt.annotate(
#             materials[i][0],
#             [x_i,ymin+0.5],
#             fontsize=20,
#             ha='center'
#         )
        
#         n=n+no_of_defects+1
    
#     if fermi_line is not None:
#         plt.axhline(y=fermi_line, color='tab:cyan', linestyle='--', alpha=0.5)

#     n=1
#     for i in range(len(materials)):
#         for row in defects:
#             if row[1] == materials[i][0]:

#                 used_y = []
#                 for j in range(len(row[2])):
#                     y_value = materials[i][1]+row[2][j][1]
#                     offset = 0
#                     for used in used_y:
#                         if abs(y_value-used) < 0.1:
#                             offset += 0.1
                    
#                     plt.hlines(y_value, n-0.1, n+0.1, lw=3)
#                     plt.annotate(
#                         f'{row[2][j][0]}',
#                         [n,y_value+offset],
#                         va='center',
#                         xytext=(30,0),
#                         textcoords='offset points',
#                         fontsize=20
#                     )
#                     used_y.append(y_value+offset)
#                 n=n+1
#         n=n+1
            
#     plt.savefig('Multi_Mat_CTLs.png')
#     plt.show()









# def plot_multi_mat_ctl_errors(defects,materials,fermi_line=None):
#     """
#     Plots all given CTLs in given materials with band alignment.

#     nb
#         Ensure that the VBM and CBM energies are given relative to one another between different materials.
#         Ensure that the 'Material Name' is the same in the defects and materials lists.
#         CTL level energies to be given relative to the VBM of the material.

#     args:

#     defects:
#             Each row: [Defect Name, Material Name, List of CTL tuplesenergies [Name, Ave_E, Min_E, Max_E] ]
#             e.g.
#                 species = [
#                     [f'N$_{'i'}$', f'HfO$_2$', [
#                         ['+1/0', 1.3469378771661913],
#                         ['0/-1', 2.8837208880316783],
#                         ['-1/-2', 2.369073993248366],
#                         ['-2/-3', 2.6696762876888496],
#                     ]],
#                     [f'N$_{'O'}$', f'HfO$_2$', [
#                         ['3+/2+', 1.7950036300732979],
#                         ['2+/1+', 3.0614261580351125],
#                         ['+1/0', 2.0040361516655647],
#                         ['0/-1', 2.164946988290245],
#                     ]],
#                     [f'N$_{'2'}$$_{'i'}$', f'HfO$_2$', [
#                         ['0/-1', 3.380293837285865],
#                         ['-1/-2', 4.198164515355801],
#                         ['-2/-3', 4.441756959205842],
#                         ['-3/-4', 3.7187798290658316]
#                     ]]
#                 ]
#     materials:
#                 Each row: [Material Name, VBM, CBM]
#                 e.g.
#                     materials = [
#                         ['Si', Sivbm, Sicbm],
#                         [f'HfO$_{'2'}$', HfO2vbm, HfO2cbm]
#                     ]
#     fermi_line:
#                 Draws dashed line at desired fermi level
#     """

#     #color_list = ["tab:orange", "tab:pink", "tab:blue", "tab:brown", "tab:red", "tab:green"]

#     color_list = {
#         '3+/2+': 'tab:red',
#         '2+/1+': 'tab:olive',
#         '1+/0': 'tab:cyan',
#         '0/1-': 'tab:blue',
#         '1-/2-': 'tab:pink',
#         '2-/3-': 'tab:purple',
#         '3-/4-': 'tab:brown'
#     }

#     defects_list = []
#     materials_list = []

#     for row in defects:
#         defects_list.append(row[:2])
#     for i in range(len(materials)):
#         materials_list.append(materials[i][0])

#     ymin = round(min(materials, key=lambda x: x[1])[1]-1)
#     ymax = round(max(materials, key=lambda x: x[2])[2]+1)

#     xmax = len(defects)+len(materials)

#     plt.figure(figsize=(30,20))
#     plt.title('Charge Transition Levels', size=60, pad=60)

#     plt.ylim(ymin,ymax)
#     plt.xlim(0,xmax)

#     n=0
#     for i in range(len(materials)):
#         no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
#         defects_list.insert(n+no_of_defects,['','',''])
#         n=n+no_of_defects+1
#     defects_list.insert(0,['','',''])

#     plt.xticks(np.linspace(0,xmax,xmax+1), [row[0] for row in defects_list], fontsize=30)
#     plt.yticks(fontsize=30)

#     plt.ylabel('Fermi Energy (eV)', size=40, labelpad=20)
#     plt.xlabel('Defects', size=40, labelpad=30) 

#     n=0
#     for i in range(len(materials)):

#         no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
#         xmin_i = n
#         xmax_i = n+no_of_defects+1
#         vbm_i = materials[i][1]
#         cbm_i = materials[i][2]

#         # plt.hlines(
#         #     y=vbm_i,
#         #     xmin=xmin_i,
#         #     xmax=xmax_i,
#         #     color='tab:green',
#         #     linestyle='-',
#         #     alpha=0.5
#         # )
        
#         # plt.hlines(
#         #     y=cbm_i,
#         #     xmin=xmin_i,
#         #     xmax=xmax_i,
#         #     color='tab:orange',
#         #     linestyle='-',
#         #     alpha=0.5
#         # )

#         plt.fill(
#             [xmin_i,xmin_i,xmax_i,xmax_i],
#             [ymin,vbm_i,vbm_i,ymin],
#             'tab:green',
#             [xmin_i,xmin_i,xmax_i,xmax_i],
#             [cbm_i,ymax,ymax,cbm_i],
#             'tab:orange',
#             alpha=0.2
#         )

#         x_i = n + (no_of_defects+1)/2

#         plt.annotate(
#             materials[i][0],
#             [x_i,ymin+0.5],
#             fontsize=30,
#             ha='center'
#         )
        
#         n=n+no_of_defects+1
    
#     if fermi_line is not None:
#         plt.axhline(y=fermi_line, color='tab:cyan', linestyle='--', alpha=0.5)

#     n=1
#     for i in range(len(materials)):
#         for row in defects:
#             if row[1] == materials[i][0]:

#                 used_y = []

#                 for j in range(len(row[2])):
#                     w = 0.4/(len(row[2])-1)
#                     m = n + w * (j - (len(row[2]) - 1) / 2)

#                     min_y = materials[i][1]+row[2][j][2]
#                     max_y = materials[i][1]+row[2][j][3]
#                     plt.fill(
#                         [m-0.05,m+0.05,m+0.05,m-0.05],
#                         [min_y,min_y, max_y,max_y],
#                         color_list.get(f'{row[2][j][0]}'),
#                         alpha=0.2
#                     )

#                 for j in range(len(row[2])):
#                     w = 0.4/(len(row[2])-1)
#                     m = n + w * (j - (len(row[2]) - 1) / 2)

#                     min_y = materials[i][1]+row[2][j][2]
#                     max_y = materials[i][1]+row[2][j][3]
#                     y_value = materials[i][1]+row[2][j][1]
#                     offset = 0
#                     for used in used_y:
#                         if abs(y_value-used) < 0.1:
#                             offset += 0.1
                    
#                     o = n-0.4 if m<n else n+0.4 if m>n else n+0.4
#                     plt.hlines(min_y, m-0.075, m+0.075, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
#                     plt.hlines(max_y, m-0.075, m+0.075, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
#                     plt.hlines(y_value, m-0.05, m+0.05, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
#                     # plt.annotate(
#                     #     f'{row[2][j][0]}',
#                     #     [o,y_value+offset],
#                     #     va='center',
#                     #     ha='center',
#                     #     xytext=(0,0),
#                     #     textcoords='offset points',
#                     #     fontsize=20,
#                     #     color=color_list.get(f'{row[2][j][0]}')
#                     # )
#                     used_y.append(y_value+offset)
#                 n=n+1
#         n=n+1
    
#     custom_leg = [[],[]]

#     for key in color_list:
#         custom_leg[0].append(key)
#         custom_leg[1].append(Line2D([0], [0], color=color_list[key], lw=6))
    
    
#     plt.legend(custom_leg[1], custom_leg[0], loc='center right', fontsize=30)

#     plt.savefig('Multi_Mat_CTLs_error.png', bbox_inches='tight')
#     plt.show()






# def temp_plot_multi_mat_ctl_errors(defects,materials,fermis):
#     """
#     Plots all given CTLs in given materials with band alignment.

#     nb
#         Ensure that the VBM and CBM energies are given relative to one another between different materials.
#         Ensure that the 'Material Name' is the same in the defects and materials lists.
#         CTL level energies to be given relative to the VBM of the material.

#     args:

#     defects:
#             Each row: [Defect Name, Material Name, List of CTL tuplesenergies [Name, Ave_E, Min_E, Max_E] ]
#             e.g.
#                 species = [
#                     [f'N$_{'i'}$', f'HfO$_2$', [
#                         ['+1/0', 1.3469378771661913],
#                         ['0/-1', 2.8837208880316783],
#                         ['-1/-2', 2.369073993248366],
#                         ['-2/-3', 2.6696762876888496],
#                     ]],
#                     [f'N$_{'O'}$', f'HfO$_2$', [
#                         ['3+/2+', 1.7950036300732979],
#                         ['2+/1+', 3.0614261580351125],
#                         ['+1/0', 2.0040361516655647],
#                         ['0/-1', 2.164946988290245],
#                     ]],
#                     [f'N$_{'2'}$$_{'i'}$', f'HfO$_2$', [
#                         ['0/-1', 3.380293837285865],
#                         ['-1/-2', 4.198164515355801],
#                         ['-2/-3', 4.441756959205842],
#                         ['-3/-4', 3.7187798290658316]
#                     ]]
#                 ]
#     materials:
#                 Each row: [Material Name, VBM, CBM]
#                 e.g.
#                     materials = [
#                         ['Si', Sivbm, Sicbm],
#                         [f'HfO$_{'2'}$', HfO2vbm, HfO2cbm]
#                     ]
#     fermi_line:
#                 Draws dashed line at desired fermi level
#     """

#     #color_list = ["tab:orange", "tab:pink", "tab:blue", "tab:brown", "tab:red", "tab:green"]

#     color_list = {
#         #'3+/2+': 'tab:red',
#         '2+/1+': 'tab:olive',
#         '1+/0': 'tab:cyan',
#         '0/1-': 'tab:blue',
#         '1-/2-': 'tab:pink',
#         #'2-/3-': 'tab:purple',
#         #'3-/4-': 'tab:brown'
#     }

#     defects_list = []
#     materials_list = []

#     for row in defects:
#         defects_list.append(row[:2])
#     for i in range(len(materials)):
#         materials_list.append(materials[i][0])

#     ymin = round(min(materials, key=lambda x: x[1])[1]-1)
#     ymax = round(max(materials, key=lambda x: x[2])[2]+1)

#     xmax = len(defects)+len(materials)

#     plt.figure(figsize=(20,20))
#     #plt.title('Charge Transition Levels', size=60, pad=60)

#     plt.ylim(ymin,ymax)
#     plt.xlim(0,xmax)

#     n=0
#     for i in range(len(materials)):
#         no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
#         defects_list.insert(n+no_of_defects,['','',''])
#         n=n+no_of_defects+1
#     defects_list.insert(0,['','',''])

#     plt.xticks(np.linspace(0,xmax,xmax+1), [row[0] for row in defects_list], fontsize=30)
#     plt.yticks(fontsize=30)

#     plt.ylabel('Fermi Energy (eV)', size=40, labelpad=20)
#     #plt.xlabel('Defects', size=40, labelpad=30) 

#     n=0
#     for i in range(len(materials)):

#         no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
#         xmin_i = n
#         xmax_i = n+no_of_defects+1
#         vbm_i = materials[i][1]
#         cbm_i = materials[i][2]

#         # plt.hlines(
#         #     y=vbm_i,
#         #     xmin=xmin_i,
#         #     xmax=xmax_i,
#         #     color='tab:green',
#         #     linestyle='-',
#         #     alpha=0.5
#         # )
        
#         # plt.hlines(
#         #     y=cbm_i,
#         #     xmin=xmin_i,
#         #     xmax=xmax_i,
#         #     color='tab:orange',
#         #     linestyle='-',
#         #     alpha=0.5
#         # )

#         plt.fill(
#             [xmin_i,xmin_i,xmax_i,xmax_i],
#             [ymin,vbm_i,vbm_i,ymin],
#             'tab:green',
#             [xmin_i,xmin_i,xmax_i,xmax_i],
#             [cbm_i,ymax,ymax,cbm_i],
#             'tab:orange',
#             alpha=0.2
#         )

#         x_i = n + (no_of_defects+1)/2

#         plt.annotate(
#             materials[i][0],
#             [x_i,ymin+0.5],
#             fontsize=30,
#             ha='center'
#         )
        
#         n=n+no_of_defects+1
    
#     for n in range(len(fermis)):
#         plt.axhline(y=fermis[n][1], label=fermis[n][0], linestyle='--', lw=3, alpha=0.8)

#     n=1
#     for i in range(len(materials)):
#         for row in defects:
#             if row[1] == materials[i][0]:

#                 used_y = []

#                 for j in range(len(row[2])):
#                     w = 0.4/(len(row[2])-1)
#                     m = n + w * (j - (len(row[2]) - 1) / 2)

#                     min_y = materials[i][1]+row[2][j][2]
#                     max_y = materials[i][1]+row[2][j][3]
#                     plt.fill(
#                         [m-0.05,m+0.05,m+0.05,m-0.05],
#                         [min_y,min_y, max_y,max_y],
#                         color_list.get(f'{row[2][j][0]}'),
#                         alpha=0.2
#                     )

#                 for j in range(len(row[2])):
#                     w = 0.4/(len(row[2])-1)
#                     m = n + w * (j - (len(row[2]) - 1) / 2)

#                     min_y = materials[i][1]+row[2][j][2]
#                     max_y = materials[i][1]+row[2][j][3]
#                     y_value = materials[i][1]+row[2][j][1]
#                     offset = 0
#                     for used in used_y:
#                         if abs(y_value-used) < 0.1:
#                             offset += 0.1
                    
#                     o = n-0.4 if m<n else n+0.4 if m>n else n+0.4
#                     plt.hlines(min_y, m-0.075, m+0.075, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
#                     plt.hlines(max_y, m-0.075, m+0.075, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
#                     plt.hlines(y_value, m-0.05, m+0.05, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
#                     # plt.annotate(
#                     #     f'{row[2][j][0]}',
#                     #     [o,y_value+offset],
#                     #     va='center',
#                     #     ha='center',
#                     #     xytext=(0,0),
#                     #     textcoords='offset points',
#                     #     fontsize=20,
#                     #     color=color_list.get(f'{row[2][j][0]}')
#                     # )
#                     used_y.append(y_value+offset)
#                 n=n+1
#         n=n+1
    
#     # custom_leg = [[],[]]

#     # for key in color_list:
#     #     custom_leg[0].append(key)
#     #     custom_leg[1].append(Line2D([0], [0], color=color_list[key], lw=6))
    
    
#     #plt.legend(custom_leg[1], custom_leg[0], loc='center right', fontsize=30)
#     #plt.legend(fontsize=20)
#     plt.savefig('Mats_fermis.png', bbox_inches='tight')
#     plt.show()