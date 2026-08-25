"""Chromatogram figures in the colours of the talk (task R3 of round 2).

One panel per fractionation window, written as a separate vector PDF so that two
panels can be placed side by side on a slide with identical axes. Each panel
shows

* the nominal densities of all species,
* the envelopes,
* the fractionation window as a shaded band,
* the worst-case densities of all species, obtained as the worst case of least
  total variation among all optimal ones (``inner_lp.smoothest_worst_case``).

All curves are drawn in pA, i.e. densities are multiplied by the mass of the
respective peak, so that they are comparable with the measured signal.

Colours are the FAU Nat palette of the talk: green is the decision, red the
species to be collected, grey what is not wanted, blue the adversary.

Usage::

    python3.12 slide_plot.py --window nominal:3.1198:3.4760 --out graphics/nominal.pdf
    python3.12 slide_plot.py --window robust:3.2009:3.3813  --out graphics/robust.pdf
"""

import argparse
import json
import os

import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import inner_lp
import run_funktionen as rf
from instanz import baue_instanz, Q_FACTOR
from params import Params

# FAU Nat palette, as used in the beamer theme of the talk
NAT_GREEN = '#43B02A'
NAT_DARK_GREEN = '#228848'
RED = '#E4322B'
GREY = '#737373'
DARK = '#3C3C3C'
BLUE = '#002F6C'

WINDOW_ALPHA = 0.14
WORST_CASE_ALPHA = 0.55


def _style(font_pt=9.5):
    """Font sizes for a panel shown at about 14 cm width."""
    plt.rcParams.update({
        'font.size': font_pt,
        'axes.labelsize': font_pt,
        'axes.titlesize': font_pt,
        'xtick.labelsize': font_pt - 0.5,
        'ytick.labelsize': font_pt - 0.5,
        'legend.fontsize': font_pt - 1.5,
        'axes.edgecolor': DARK,
        'axes.labelcolor': DARK,
        'xtick.color': DARK,
        'ytick.color': DARK,
        'text.color': DARK,
        'pdf.fonttype': 42,
        'savefig.transparent': True,
    })


def compute(inst, lo, hi, tol=1e-9, objective='tv', mip_indicator=True):
    """Worst-case measures of least total variation for all species."""
    return [inner_lp.smoothest_worst_case(inst, i, lo, hi, tol=tol, objective=objective,
                                          mip_indicator=mip_indicator)
            for i in inst.groessen]


def panel(inst, measures, lo, hi, path, desired, xlim=None, ylim=None,
          show_contaminant_envelopes=True, width_cm=14., height_cm=8.4, font_pt=9.5,
          smoothed=True, legend=True):
    """Write one panel."""
    _style(font_pt)
    t = np.asarray(inst.time_points)
    K = inst.anzahl_prozess
    d = inst.zeit_diskret
    # signal in pA: the normalised densities are scaled back by the peak mass
    scale = [inst.q0[i] / Q_FACTOR for i in inst.groessen]

    fig, ax = plt.subplots(figsize=(width_cm / 2.54, height_cm / 2.54))

    # fractionation window
    ax.axvspan(t[lo], t[hi], facecolor=NAT_GREEN, alpha=WINDOW_ALPHA, lw=0., zorder=0)
    for x in (t[lo], t[hi]):
        ax.axvline(x, color=NAT_DARK_GREEN, lw=1.0, zorder=1)

    key = 'cell_mass_smoothed' if smoothed else 'cell_mass_vertex'
    for i in inst.groessen:
        is_desired = (i == desired)
        nominal = np.asarray(inst.matrix_nom[i]) * scale[i]
        envelope = np.asarray(inst.schlauch_rtd[i]) * scale[i]
        worst = np.asarray(measures[i][key]) / d * scale[i]

        # worst case, filled to zero
        ax.fill_between(t[:K], 0., worst, step='post',
                        color=RED if is_desired else BLUE,
                        alpha=WORST_CASE_ALPHA, lw=0., zorder=2)
        ax.plot(t[:K], worst, drawstyle='steps-post',
                color=RED if is_desired else BLUE, lw=0.7,
                alpha=WORST_CASE_ALPHA, zorder=3)
        # envelope
        if is_desired or show_contaminant_envelopes:
            ax.plot(t, envelope, color=RED if is_desired else DARK, lw=0.8,
                    ls=(0, (4, 2)), zorder=4)
        # nominal density
        ax.plot(t, nominal, color=RED if is_desired else GREY,
                lw=1.6 if is_desired else 0.9, zorder=5)

    ax.set_xlabel('Time in min')
    ax.set_ylabel('Signal in pA')
    ax.set_xlim(*(xlim if xlim else (t[0], t[-1])))
    if ylim:
        ax.set_ylim(*ylim)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    if legend:
        _legend(ax, show_contaminant_envelopes)

    fig.tight_layout(pad=0.3)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    print(f"written: {path}")


def _legend(ax, show_contaminant_envelopes):
    handles = [
            Patch(facecolor=NAT_GREEN, alpha=WINDOW_ALPHA, edgecolor=NAT_DARK_GREEN,
                  label='fractionation window'),
            Line2D([0], [0], color=RED, lw=1.6, label='desired, nominal'),
            Line2D([0], [0], color=RED, lw=0.8, ls=(0, (4, 2)), label='desired, envelope'),
            Patch(facecolor=RED, alpha=WORST_CASE_ALPHA, label='desired, worst case'),
        Line2D([0], [0], color=GREY, lw=0.9, label='contaminants, nominal'),
    ]
    if show_contaminant_envelopes:
        handles.append(Line2D([0], [0], color=DARK, lw=0.8, ls=(0, (4, 2)),
                              label='contaminants, envelope'))
    handles.append(Patch(facecolor=BLUE, alpha=WORST_CASE_ALPHA,
                         label='contaminants, worst case'))
    ax.legend(handles=handles, loc='upper right', frameon=True, framealpha=0.95,
              edgecolor=DARK, fancybox=False, borderpad=0.5, handlelength=1.6,
              labelspacing=0.35).get_frame().set_linewidth(0.6)


def panel_from_json(payload, path, xlim=None, ylim=None, desired=2,
                    show_contaminant_envelopes=True, width_cm=14., height_cm=8.4,
                    font_pt=9.5, smoothed=True, legend=True):
    """Re-draw a panel from the JSON written by a previous run, without Gurobi.

    Useful for changing sizes, fonts or limits without re-solving the linear
    programs; all curves are stored in the JSON in pA.
    """
    _style(font_pt)
    t = np.asarray(payload['time_points'])
    d = payload['delta_N']
    lo, hi = payload['index_lower'], payload['index_upper']
    key = 'worst_case_pA' if smoothed else 'worst_case_vertex_pA'

    fig, ax = plt.subplots(figsize=(width_cm / 2.54, height_cm / 2.54))
    ax.axvspan(t[lo], t[hi], facecolor=NAT_GREEN, alpha=WINDOW_ALPHA, lw=0., zorder=0)
    for x in (t[lo], t[hi]):
        ax.axvline(x, color=NAT_DARK_GREEN, lw=1.0, zorder=1)

    for i, (nominal, envelope, worst) in enumerate(zip(payload['nominal_pA'],
                                                       payload['envelope_pA'],
                                                       payload[key])):
        is_desired = (i == desired)
        worst = np.asarray(worst)
        ax.fill_between(t[:len(worst)], 0., worst, step='post',
                        color=RED if is_desired else BLUE, alpha=WORST_CASE_ALPHA,
                        lw=0., zorder=2)
        ax.plot(t[:len(worst)], worst, drawstyle='steps-post',
                color=RED if is_desired else BLUE, lw=0.7, alpha=WORST_CASE_ALPHA, zorder=3)
        if is_desired or show_contaminant_envelopes:
            ax.plot(t, envelope, color=RED if is_desired else DARK, lw=0.8,
                    ls=(0, (4, 2)), zorder=4)
        ax.plot(t, nominal, color=RED if is_desired else GREY,
                lw=1.6 if is_desired else 0.9, zorder=5)

    ax.set_xlabel('Time in min')
    ax.set_ylabel('Signal in pA')
    ax.set_xlim(*(xlim if xlim else (t[0], t[-1])))
    if ylim:
        ax.set_ylim(*ylim)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    if legend:
        _legend(ax, show_contaminant_envelopes)
    fig.tight_layout(pad=0.3)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    print(f"written: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--sample', default='long')
    parser.add_argument('--aggregation_factor', type=int, default=1)
    parser.add_argument('--reinheit', type=float, default=.95)
    parser.add_argument('--wunschgroesse', type=int, default=2)
    parser.add_argument('--from_json', default=None,
                        help="Re-draw from the JSON of a previous run instead of solving.")
    parser.add_argument('--window', default=None,
                        help="name:x_lower:x_upper, e.g. robust:3.2009:3.3813")
    parser.add_argument('--out', required=True)
    parser.add_argument('--json', default=None)
    parser.add_argument('--tol', type=float, default=1e-9)
    parser.add_argument('--objective', default='tv', choices=('tv', 'curvature'))
    parser.add_argument('--vertex', action='store_true',
                        help="Plot the unsmoothed vertex solution instead.")
    parser.add_argument('--no_contaminant_envelopes', action='store_true')
    parser.add_argument('--xlim', type=float, nargs=2, default=[2.85, 3.62])
    parser.add_argument('--ylim', type=float, nargs=2, default=[0., 12.5])
    parser.add_argument('--width_cm', type=float, default=14.)
    parser.add_argument('--height_cm', type=float, default=8.4)
    parser.add_argument('--font_pt', type=float, default=9.5)
    parser.add_argument('--no_legend', action='store_true')
    args = parser.parse_args()

    if args.from_json:
        with open(args.from_json) as f:
            payload = json.load(f)
        print(f"{payload['window_name']}: replotting from {args.from_json}")
        panel_from_json(payload, args.out, xlim=args.xlim, ylim=args.ylim,
                        desired=args.wunschgroesse,
                        show_contaminant_envelopes=not args.no_contaminant_envelopes,
                        width_cm=args.width_cm, height_cm=args.height_cm,
                        font_pt=args.font_pt, smoothed=not args.vertex,
                        legend=not args.no_legend)
        return

    if not args.window:
        parser.error("either --window or --from_json is required")
    name, xl, xu = args.window.split(':')
    params = Params(aggregation_factor=args.aggregation_factor, reinheit=args.reinheit,
                    wunschgroesse=args.wunschgroesse, sample=args.sample, fix=False,
                    fix_lower=0, fix_upper=0, nominal=False, plot=False)
    data = rf.read_data(params)
    inst = baue_instanz(data[0], data[1], data[2], data[3], params)
    t = np.asarray(inst.time_points)
    lo = int(round((float(xl) - t[0]) / inst.zeit_diskret))
    hi = int(round((float(xu) - t[0]) / inst.zeit_diskret))
    print(f"{name}: window [{t[lo]:.5f}, {t[hi]:.5f}], indices ({lo}, {hi}), "
          f"delta_N = {inst.zeit_diskret:.6g}")

    measures = compute(inst, lo, hi, tol=args.tol, objective=args.objective)
    for m in measures:
        print(f"  species index {m['species_index']}: v* = {m['v_star']:.9f}, smoothed "
              f"{m['v_smoothed']:.9f} (deviation {m['value_deviation']:.2e}), "
              f"P(window) = {m['prob_in_window']:.6f}, total variation "
              f"{m['total_variation_vertex']:.1f} -> {m['total_variation_smoothed']:.1f} "
              f"({100 * m['total_variation_reduction']:.1f} % less), cells with mass "
              f"{m['n_cells_with_mass']}, of those at the envelope "
              f"{m['n_cells_saturating_envelope']} "
              f"({100 * m['fraction_saturating']:.1f} %)")

    panel(inst, measures, lo, hi, args.out, desired=args.wunschgroesse,
          xlim=args.xlim, ylim=args.ylim,
          show_contaminant_envelopes=not args.no_contaminant_envelopes,
          width_cm=args.width_cm, height_cm=args.height_cm, font_pt=args.font_pt,
          smoothed=not args.vertex, legend=not args.no_legend)

    if args.json:
        payload = {
            'window_name': name, 'x_lower': float(t[lo]), 'x_upper': float(t[hi]),
            'index_lower': lo, 'index_upper': hi, 'delta_N': inst.zeit_diskret,
            'tolerance': args.tol, 'smoothing_objective': args.objective,
            'plotted': 'vertex' if args.vertex else 'smoothed',
            'colours': {'window_fill': NAT_GREEN, 'window_edge': NAT_DARK_GREEN,
                        'desired': RED, 'contaminant_nominal': GREY,
                        'contaminant_envelope': DARK, 'contaminant_worst_case': BLUE,
                        'axes': DARK},
            'q0_scaled': list(inst.q0), 'q_factor': Q_FACTOR,
            'species': [{k: v for k, v in m.items()
                         if k not in ('cell_mass_vertex', 'cell_mass_smoothed', 'cap')}
                        for m in measures],
            'time_points': [float(v) for v in t],
            'nominal_pA': [[float(v) for v in np.asarray(inst.matrix_nom[i])
                            * inst.q0[i] / Q_FACTOR] for i in inst.groessen],
            'envelope_pA': [[float(v) for v in np.asarray(inst.schlauch_rtd[i])
                             * inst.q0[i] / Q_FACTOR] for i in inst.groessen],
            'worst_case_pA': [[float(v) for v in np.asarray(m['cell_mass_smoothed'])
                               / inst.zeit_diskret * inst.q0[i] / Q_FACTOR]
                              for i, m in zip(inst.groessen, measures)],
            'worst_case_vertex_pA': [[float(v) for v in np.asarray(m['cell_mass_vertex'])
                                      / inst.zeit_diskret * inst.q0[i] / Q_FACTOR]
                                     for i, m in zip(inst.groessen, measures)],
        }
        with open(args.json, 'w') as f:
            json.dump(payload, f, default=float)
        print(f"written: {args.json}")


if __name__ == '__main__':
    main()
