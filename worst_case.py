"""Worst-case measures of the safe approximation (Remark "no_atoms").

For a given fractionation interval the adversarial problem of every species is
solved as the linear program "safe_approx_primal", see ``inner_lp.LowerLP``. Its
variables are exactly the weights ``w^-_tau`` and ``w^+_tau`` of the re-dualised
problem, i.e. the dual variables of Constraints (32c) and (32d) of the
mixed-integer model. They are read off the primal LP directly instead of from the
duals of the MIP, which is equivalent (the LP is the re-dualisation of the same
system) and avoids fixing the binaries of the MIP and re-solving it.

The script writes a figure with one panel per selected species showing

* the nominal density,
* the envelope,
* the density ``(w^-_tau + w^+_tau) / delta_N`` of the worst-case grid measure,

together with the fractionation window. The last quantity is the object referred
to in Remark "no_atoms": the mass of a single grid point is bounded by the
envelope mass of a cell of length ``delta_N``, so that the plotted density stays
below the envelope and no atoms appear in the limit.

Usage::

    python3.12 worst_case.py --x_lower 3.2009 --x_upper 3.3813 \\
        --out Plots/worst_case_measures.pdf
"""

import argparse
import json
import os

import numpy as np
from matplotlib import pyplot as plt

import inner_lp
import run_funktionen as rf
from instanz import baue_instanz
from params import Params


def worst_case_measures(inst, lo, hi, mip_indicator=True):
    """One worst-case measure per species for the interval ``[t[lo], t[hi]]``."""
    lps = inner_lp.build_all(inst, 'lower', mip_indicator=mip_indicator)
    out = []
    for i, lp in enumerate(lps):
        val = lp.solve(lo, hi)
        K = lp.K
        wm = np.array([lp.wm[j].X for j in range(K)])
        wp = np.zeros(K)
        wp[:K - 1] = np.array([lp.wp[j].X for j in range(K - 1)])
        cap = inst.zeit_diskret * np.asarray(inst.schlauch_rtd[i][:K])
        out.append({'species_index': i, 'val': val, 'wm': wm, 'wp': wp,
                    'cell_mass': wm + wp, 'cap': cap})
    return out


def _binned(values, k):
    """Sum of ``values`` over blocks of ``k`` entries, padded to full blocks."""
    n = int(np.ceil(len(values) / k)) * k
    padded = np.zeros(n)
    padded[:len(values)] = values
    return padded.reshape(-1, k).sum(axis=1)


def plot(inst, measures, lo, hi, species, path, labels=None, title=None, bin_cells=50,
         xlim=None):
    """Figure with one panel per species.

    The optimal solution of the adversarial LP is highly degenerate: every
    selection of cells that saturates the envelope and meets the moment
    constraints has the same objective value, so the cell-by-cell structure of
    the worst-case measure is arbitrary. It is therefore drawn twice, once
    cell by cell (thin) and once averaged over ``bin_cells`` cells (thick). The
    averaged curve is the density of the measure P_1 of the convergence proof,
    i.e. the envelope scaled by the mass that the worst case places in the
    respective part of the domain.
    """
    t = np.asarray(inst.time_points)
    K = inst.anzahl_prozess
    d = inst.zeit_diskret
    fig, axes = plt.subplots(len(species), 1, figsize=(7.2, 2.7 * len(species)), sharex=True)
    if len(species) == 1:
        axes = [axes]
    for ax, i in zip(axes, species):
        mrec = measures[i]
        env = np.asarray(inst.schlauch_rtd[i])
        ax.plot(t, inst.matrix_nom[i], color='0.55', lw=1.0, label='nominal density')
        ax.plot(t, env, color='black', lw=1.0, label='envelope')
        ax.plot(t[:K], mrec['cell_mass'] / d, color='firebrick', lw=0.4, alpha=0.35,
                label=r'worst case, cell by cell: $(w^-_\tau+w^+_\tau)/\delta_N$')
        binned = _binned(mrec['cell_mass'], bin_cells) / (bin_cells * d)
        edges = t[0] + d * bin_cells * np.arange(len(binned) + 1)
        ax.stairs(binned, edges, color='firebrick', lw=1.4,
                  label=f'worst case, averaged over {bin_cells} cells')
        ax.axvspan(t[lo], t[hi], color='tab:blue', alpha=0.18, lw=0,
                   label='fractionation window')
        name = labels[i] if labels else f'species index {i}'
        ax.set_title(f"{name}: worst case collects "
                     f"{mrec['val'] / inst.a_coeff(i) * 100:.1f} % of its mass in the window",
                     fontsize=9)
        ax.set_ylabel('density (1/min)', fontsize=8)
        ax.tick_params(labelsize=8)
        ax.set_xlim(*(xlim if xlim else (t[0], t[-1])))
    axes[-1].set_xlabel('time in min', fontsize=8)
    axes[0].legend(fontsize=7, loc='upper left', framealpha=0.9)
    if title:
        fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, bbox_inches='tight')
    print(f"written: {path}")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sample', default='long')
    parser.add_argument('--aggregation_factor', type=int, default=1)
    parser.add_argument('--reinheit', type=float, default=.95)
    parser.add_argument('--wunschgroesse', type=int, default=2)
    parser.add_argument('--x_lower', type=float, default=3.2009)
    parser.add_argument('--x_upper', type=float, default=3.3813)
    parser.add_argument('--species', type=int, nargs='+', default=[2, 3],
                        help="Indices of the species to plot. Default 2 (desired, s=32) and "
                             "3 (contaminant s=33).")
    parser.add_argument('--names', nargs='+', default=['s=30', 's=31', 's=32 (desired)', 's=33'])
    parser.add_argument('--bin_cells', type=int, default=50,
                        help="Number of cells the worst-case density is additionally averaged "
                             "over. Default 50, i.e. 0.005 min.")
    parser.add_argument('--xlim', type=float, nargs=2, default=None,
                        help="Restrict the plotted time interval.")
    parser.add_argument('--out', default='output/worst_case_measures.pdf')
    parser.add_argument('--json', default=None)
    parser.add_argument('--theory_indicator', action='store_true',
                        help="Use the indicator 1^c of the convergence proof instead of the "
                             "weakened indicator of the mixed-integer model.")
    args = parser.parse_args()

    params = Params(aggregation_factor=args.aggregation_factor, reinheit=args.reinheit,
                    wunschgroesse=args.wunschgroesse, sample=args.sample, fix=False,
                    fix_lower=0, fix_upper=0, nominal=False, plot=False)
    data = rf.read_data(params)
    inst = baue_instanz(data[0], data[1], data[2], data[3], params)
    t = np.asarray(inst.time_points)
    lo = int(round((args.x_lower - t[0]) / inst.zeit_diskret))
    hi = int(round((args.x_upper - t[0]) / inst.zeit_diskret))
    print(f"grid: {len(t)} points, delta_N = {inst.zeit_diskret:.6g}; "
          f"interval [{t[lo]:.5f}, {t[hi]:.5f}]")

    measures = worst_case_measures(inst, lo, hi, mip_indicator=not args.theory_indicator)

    summary = []
    for i in inst.groessen:
        mrec = measures[i]
        cm, cap = mrec['cell_mass'], mrec['cap']
        active = cm > 1e-12
        saturated = cm > 0.999 * cap
        rec = {
            'species_index': i,
            'val': float(mrec['val']),
            'prob_in_window': float(mrec['val'] / inst.a_coeff(i)),
            'n_cells_with_mass': int(active.sum()),
            'n_cells_saturating_envelope': int((saturated & active).sum()),
            'max_cell_mass': float(cm.max()),
            'max_envelope_cell_mass': float(cap.max()),
            'max_ratio_mass_to_cap': float((cm[active] / cap[active]).max()) if active.any() else 0.,
            'mean': float(np.sum(mrec['wm'] * t[:inst.anzahl_prozess]
                                 + mrec['wp'] * t[1:inst.anzahl_prozess + 1])),
        }
        summary.append(rec)
        print(f"  species index {i}: val={rec['val']:.6f}, P(window)={rec['prob_in_window']:.6f}, "
              f"cells with mass {rec['n_cells_with_mass']}, saturating the envelope "
              f"{rec['n_cells_saturating_envelope']}, max cell mass {rec['max_cell_mass']:.3e} "
              f"(<= {rec['max_envelope_cell_mass']:.3e}), mean {rec['mean']:.6f}")

    plot(inst, measures, lo, hi, args.species, args.out, labels=args.names,
         bin_cells=args.bin_cells, xlim=args.xlim)
    if args.json:
        with open(args.json, 'w') as f:
            json.dump({'delta_N': inst.zeit_diskret, 'x_lower': float(t[lo]),
                       'x_upper': float(t[hi]), 'mip_indicator': not args.theory_indicator,
                       'species': summary}, f, indent=2)
        print(f"written: {args.json}")


if __name__ == "__main__":
    main()
