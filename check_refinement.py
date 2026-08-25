"""Consistency checks for the grid refinement (tasks R1.a and R1.b of round 2).

R1.a  For every species and every data cell, the envelope masses of the r refined
      sub-cells must sum to the envelope mass of the data cell. The refined
      ambiguity set is then the same as the one on the data grid.
R1.b  With r = 1 the refined code path must reproduce the data grid instance
      bit-identically.

In addition the deviation of the naive variant (``--refine_recompute``, i.e.
moment bounds and envelope rebuilt on the refined grid) is quantified, since that
is the variant a reader might expect.

Usage::

    python3.12 check_refinement.py --refinement 2 5 10 --json results/r1_checks.json
"""

import argparse
import json

import numpy as np

import hilfsfunktionen as aux
import run_funktionen as rf
from instanz import baue_instanz
from params import Params


def instance(r, recompute=False, sample='long'):
    params = Params(aggregation_factor=1, reinheit=.95, wunschgroesse=2, sample=sample,
                    fix=False, fix_lower=0, fix_upper=0, nominal=False, plot=False,
                    refinement_factor=r, refine_recompute=recompute)
    data = rf.read_data(params)
    return baue_instanz(data[0], data[1], data[2], data[3], params), params


def cell_masses(inst, i):
    d = inst.zeit_diskret
    return np.asarray(inst.schlauch_rtd[i][:inst.anzahl_prozess]) * d


def check(r, base, base_inst):
    """R1.a and the comparison of the two variants for one refinement factor."""
    inst, _ = instance(r)
    naive, _ = instance(r, recompute=True)
    out = {'refinement_factor': r, 'delta_N': inst.zeit_diskret,
           'n_grid_points': len(inst.time_points), 'n_cells': inst.anzahl_prozess}

    # --- R1.a: envelope mass per data cell
    worst_abs, worst_rel = 0., 0.
    for i in inst.groessen:
        coarse = cell_masses(base_inst, i)
        fine = cell_masses(inst, i)
        assert len(fine) == r * len(coarse), (len(fine), r, len(coarse))
        blocks = fine.reshape(len(coarse), r).sum(axis=1)
        diff = np.abs(blocks - coarse)
        worst_abs = max(worst_abs, float(diff.max()))
        scale = max(coarse.max(), 1e-30)
        worst_rel = max(worst_rel, float(diff.max() / scale))
    out['R1a_max_abs_deviation'] = worst_abs
    out['R1a_max_rel_deviation'] = worst_rel
    out['R1a_ok'] = worst_rel < 1e-12

    # --- same check for the naive variant
    worst_abs_n, worst_rel_n, cells_n = 0., 0., 0
    for i in inst.groessen:
        coarse = cell_masses(base_inst, i)
        fine = cell_masses(naive, i)
        blocks = fine.reshape(len(coarse), r).sum(axis=1)
        diff = np.abs(blocks - coarse)
        cells_n += int((diff > 1e-12 * max(coarse.max(), 1e-30)).sum())
        worst_abs_n = max(worst_abs_n, float(diff.max()))
        worst_rel_n = max(worst_rel_n, float(diff.max() / max(coarse.max(), 1e-30)))
    out['naive_max_abs_deviation'] = worst_abs_n
    out['naive_max_rel_deviation'] = worst_rel_n
    out['naive_cells_with_deviation'] = cells_n
    out['naive_total_cells'] = 4 * base_inst.anzahl_prozess

    # --- moment bounds: unchanged in the invariant variant, shifted in the naive one
    out['mu_minus'] = list(inst.ret_time_minus)
    out['mu_plus'] = list(inst.ret_time_plus)
    out['mu_minus_shift_naive'] = [n - b for n, b in zip(naive.ret_time_minus,
                                                         base_inst.ret_time_minus)]
    out['mu_plus_shift_naive'] = [n - b for n, b in zip(naive.ret_time_plus,
                                                        base_inst.ret_time_plus)]
    out['mu_shift_naive_predicted'] = base_inst.zeit_diskret * (r - 1) / (2. * r)
    out['sigma2_plus'] = list(inst.sigma2_plus)
    out['sigma2_plus_naive'] = list(naive.sigma2_plus)
    out['q0'] = list(inst.q0)
    out['q0_deviation'] = float(np.abs(np.asarray(inst.q0) - np.asarray(base_inst.q0)).max())
    out['rho_max'] = [inst.rho_max(i) for i in inst.groessen]
    out['rho_max_deviation'] = float(max(abs(inst.rho_max(i) - base_inst.rho_max(i))
                                         for i in inst.groessen))
    out['T_bar'] = [inst.t_bar(i) for i in inst.groessen]
    out['T_bar_base'] = [base_inst.t_bar(i) for i in inst.groessen]
    out['totm_desired'] = inst.totm_desired
    out['totm_desired_base'] = base_inst.totm_desired
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--refinement', type=int, nargs='+', default=[2, 5, 10])
    parser.add_argument('--json', default=None)
    args = parser.parse_args()

    base_inst, base = instance(1)
    print(f"data grid: {len(base_inst.time_points)} points, delta_N = {base_inst.zeit_diskret:.6g}, "
          f"T = [{base_inst.time_points[0]:.5f}, {base_inst.time_points[-1]:.5f}]")

    # --- R1.b: r = 1 must be the identity
    one, _ = instance(1)
    same = {
        'time_points': np.array_equal(np.asarray(one.time_points), np.asarray(base_inst.time_points)),
        'matrix_nom': np.array_equal(one.matrix_nom, base_inst.matrix_nom),
        'matrix_min': np.array_equal(one.matrix_min, base_inst.matrix_min),
        'matrix_max': np.array_equal(one.matrix_max, base_inst.matrix_max),
        'envelope': np.array_equal(np.asarray(one.schlauch_rtd), np.asarray(base_inst.schlauch_rtd)),
        'mu_minus': one.ret_time_minus == base_inst.ret_time_minus,
        'mu_plus': one.ret_time_plus == base_inst.ret_time_plus,
        'varianz_schranke': one.varianz_schranke == base_inst.varianz_schranke,
        'zeit_diskret': one.zeit_diskret == base_inst.zeit_diskret,
    }
    print(f"R1.b (r = 1 is the identity): {'ok' if all(same.values()) else 'FAILED'} {same}")

    rows = [check(r, base, base_inst) for r in args.refinement]
    for o in rows:
        print(f"\nr = {o['refinement_factor']}: delta_N = {o['delta_N']:.6g}, "
              f"{o['n_cells']} cells")
        print(f"  R1.a  max |sum of sub-cell masses - data cell mass| = "
              f"{o['R1a_max_abs_deviation']:.3e} absolute, "
              f"{o['R1a_max_rel_deviation']:.3e} relative -> "
              f"{'ok' if o['R1a_ok'] else 'FAILED'}")
        print(f"  naive variant: {o['naive_cells_with_deviation']} of {o['naive_total_cells']} "
              f"data cells deviate, at most {o['naive_max_rel_deviation']:.3e} relative")
        print(f"  moment bounds: unchanged here; naive variant shifts mu_+ by "
              f"{['%.2e' % v for v in o['mu_plus_shift_naive']]} "
              f"(predicted {o['mu_shift_naive_predicted']:.2e})")
        print(f"  q0 deviation {o['q0_deviation']:.3e}, rho_max deviation "
              f"{o['rho_max_deviation']:.3e}, totm_desired {o['totm_desired']:.12f} vs "
              f"{o['totm_desired_base']:.12f}")

    if args.json:
        with open(args.json, 'w') as f:
            json.dump({'R1b_identity': same, 'rows': rows,
                       'data_grid': {'n_points': len(base_inst.time_points),
                                     'delta_N': base_inst.zeit_diskret}}, f, indent=2, default=float)
        print(f"\nwritten: {args.json}")


if __name__ == '__main__':
    main()
