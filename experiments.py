"""Driver for the computational experiments reported in the manuscript.

Every sub-command writes one JSON file with all numbers it produced and one
Gurobi log per solve, so that every reported number can be traced back to a raw
log. Nothing is cached: re-running a sub-command re-solves everything.

    python3.12 experiments.py c1 --out DIR   # model sizes and solver statistics
    python3.12 experiments.py c2 --out DIR   # certified enclosure (Corollary "enclosure")
    python3.12 experiments.py c4 --out DIR   # effect of the second moment constraint
    python3.12 experiments.py c5 --out DIR   # dependence on the size of the ambiguity set
    python3.12 experiments.py eval --out DIR # evaluate given fractionation intervals
    python3.12 experiments.py v1 --out DIR   # moment and variance bounds of the instances

The aggregation factors and the corresponding grid widths delta_N of the default
data set 'long' are

    factor    1     2      5      10     20     50     100
    delta_N   1e-4  2e-4   5e-4   1e-3   2e-3   5e-3   1e-2
"""

import argparse
import json
import os
import platform
import subprocess
import time

import numpy as np

import error_bound
import hilfsfunktionen as aux
import inner_lp
import run_funktionen as rf
from instanz import baue_instanz
from params import Params

# aggregation factors of Table "runtime", in the order of the table
FACTORS_TABLE = [100, 50, 20, 10, 5, 2, 1]
# the four finest discretisations, used for the enclosure
FACTORS_FINE = [10, 5, 2, 1]


def environment():
    """Hardware and software the numbers were produced on."""
    import gurobipy as gp
    import matplotlib
    cpu = None
    try:
        with open('/proc/cpuinfo') as f:
            for line in f:
                if line.startswith('model name'):
                    cpu = line.split(':', 1)[1].strip()
                    break
    except OSError:
        pass
    mem = None
    try:
        with open('/proc/meminfo') as f:
            mem = f.readline().split()[1] + ' kB'
    except OSError:
        pass
    try:
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                         stderr=subprocess.DEVNULL).decode().strip()
        dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'],
                                            stderr=subprocess.DEVNULL).decode().strip())
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    return {
        'python': platform.python_version(),
        'platform': platform.platform(),
        'cpu': cpu,
        'cpu_count': os.cpu_count(),
        'memory_total': mem,
        'gurobi': '.'.join(str(v) for v in gp.gurobi.version()),
        'numpy': np.__version__,
        'matplotlib': matplotlib.__version__,
        'git_commit': commit,
        'git_dirty': dirty,
        'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S'),
    }


def make_params(**kwargs):
    base = dict(aggregation_factor=1, reinheit=.95, wunschgroesse=2, sample='long',
                fix=False, fix_lower=0, fix_upper=0, nominal=False, plot=False)
    base.update(kwargs)
    return Params(**base)


def _logdir(out):
    logs = os.path.join(out, 'logs')
    os.makedirs(logs, exist_ok=True)
    return logs


def _dump(out, name, payload):
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, name)
    payload.setdefault('environment', environment())
    with open(path, 'w') as f:
        json.dump(payload, f, indent=2, default=float)
    print(f"\nwritten: {path}")
    return path


def _solve(out, tag, **kwargs):
    """One model solve with logging; returns the result dictionary."""
    logs = _logdir(out)
    log = os.path.join(logs, f'{tag}.log')
    if os.path.exists(log):
        os.remove(log)
    params = make_params(log_file=log, stats_file=os.path.join(logs, f'{tag}.json'), **kwargs)
    print(f"\n### {tag}: {kwargs}", flush=True)
    res = rf.run(params)
    res['log_file'] = log
    res['tag'] = tag
    return res


def _slim(res):
    """The parts of a result that go into the aggregated JSON."""
    keys = ('stats', 'delta_N', 'index_lower', 'index_upper', 'x_lower', 'x_upper',
            'interval_length', 'obj_nom', 'obj_robust', 'obj_interval', 'purity_dual',
            'purity_sampled', 'collected_mass', 'inner_values', 'a_coeff', 'log_file', 'tag',
            'running_time')
    return {k: res[k] for k in keys if k in res}


# --------------------------------------------------------------------------- C1
def task_c1(args):
    """Model sizes, running times and B&B nodes for all discretisations."""
    runs = []
    for af in FACTORS_TABLE:
        res = _solve(args.out, f'c1_af{af}', aggregation_factor=af)
        runs.append(_slim(res))
    payload = {'task': 'C1', 'sample': 'long', 'runs': runs}
    _dump(args.out, 'c1_model_sizes.json', payload)
    print(f"\n{'delta_N':>10} {'vars':>8} {'cont':>8} {'bin':>8} {'constrs':>9} "
          f"{'time(s)':>9} {'nodes':>8} {'x^-':>9} {'x^+':>9} {'length':>8}")
    for r in runs:
        s = r['stats']
        print(f"{r['delta_N']:>10.6g} {s['num_vars']:>8} {s['num_continuous_vars']:>8} "
              f"{s['num_bin_vars']:>8} {s['num_constrs']:>9} {s['runtime']:>9.2f} "
              f"{s['node_count']:>8.0f} {r.get('x_lower', float('nan')):>9.4f} "
              f"{r.get('x_upper', float('nan')):>9.4f} {r.get('interval_length', float('nan')):>8.4f}")


# --------------------------------------------------------------------------- C2
def task_c2(args):
    """Certified enclosure of Corollary "enclosure"."""
    rows = []
    for af in FACTORS_FINE:
        params = make_params(aggregation_factor=af)
        data = rf.read_data(params)
        inst = baue_instanz(data[0], data[1], data[2], data[3], params)
        bounds = error_bound.delta_bounds(inst)
        total = sum(b['Delta_N'] for b in bounds)
        print(f"\n### c2_af{af}: delta_N = {inst.zeit_diskret:.6g}, sum_s Delta_N^s = {total:.6f}")

        # The relaxed model is much harder than the original one, so it is solved
        # with a single objective (the interval length, which is the only quantity
        # the enclosure refers to) and with a time limit. If the time limit is
        # hit, the dual bound of the truncated run is still a valid upper bound on
        # val_N^up and therefore on the optimal value of the original problem.
        base = _solve(args.out, f'c2_base_af{af}', aggregation_factor=af,
                      single_objective=True)
        relaxed = _solve(args.out, f'c2_relaxed_af{af}', aggregation_factor=af,
                         purity_rhs=-total, single_objective=True,
                         time_limit=args.time_limit)

        val_n = base.get('interval_length')
        val_up = relaxed.get('interval_length')
        val_up_bound = relaxed['stats'].get('length_bound')
        row = {
            'aggregation_factor': af,
            'delta_N': inst.zeit_diskret,
            'delta_bounds': bounds,
            'total_delta': total,
            'val_N': val_n,
            'val_N_up': val_up,
            'val_N_up_proven_optimal': relaxed['stats'].get('proven_optimal'),
            'val_N_up_dual_bound': val_up_bound,
            'width': (val_up - val_n) if val_n is not None and val_up is not None else None,
            'width_certified': (val_up_bound - val_n)
            if val_n is not None and val_up_bound is not None else None,
            'base': _slim(base),
            'relaxed': _slim(relaxed),
        }
        rows.append(row)
    payload = {'task': 'C2', 'time_limit': args.time_limit, 'rows': rows}
    _dump(args.out, 'c2_enclosure.json', payload)
    print(f"\n{'delta_N':>10} {'val_N':>10} {'sum Delta':>12} {'val_N^up':>10} {'width':>10} "
          f"{'bound':>10} {'optimal':>8}")
    for r in rows:
        print(f"{r['delta_N']:>10.6g} {r['val_N']:>10.4f} {r['total_delta']:>12.4f} "
              f"{r['val_N_up']:>10.4f} {r['width']:>10.4f} "
              f"{r['val_N_up_dual_bound'] or float('nan'):>10.4f} "
              f"{str(r['val_N_up_proven_optimal']):>8}")


# --------------------------------------------------------------------------- C4
def task_c4(args):
    """Effect of the relaxed second moment constraint."""
    af = args.aggregation_factor
    with_moment = _solve(args.out, f'c4_with_af{af}', aggregation_factor=af)
    without = _solve(args.out, f'c4_without_af{af}', aggregation_factor=af, no_second_moment=True)

    params = make_params(aggregation_factor=af)
    data = rf.read_data(params)
    inst = baue_instanz(data[0], data[1], data[2], data[3], params)
    lemma = [{'species_index': i,
              'sigma2_plus': inst.sigma2_plus[i],
              'relaxation_term': ((inst.ret_time_plus[i] - inst.ret_time_minus[i]) ** 2) / 4.,
              'variance_bound_of_lemma': inst.sigma2_plus[i]
              + ((inst.ret_time_plus[i] - inst.ret_time_minus[i]) ** 2) / 4.}
             for i in inst.groessen]
    payload = {'task': 'C4', 'aggregation_factor': af, 'delta_N': inst.zeit_diskret,
               'with_second_moment': _slim(with_moment), 'without_second_moment': _slim(without),
               'lemma_bound_variance_relaxation': lemma}
    _dump(args.out, 'c4_second_moment.json', payload)
    for name, r in (('with', with_moment), ('without', without)):
        print(f"{name:>8}: interval ({r.get('x_lower')}, {r.get('x_upper')}), "
              f"length {r.get('interval_length')}, collected mass {r.get('collected_mass')}, "
              f"purity(sampled) {r.get('purity_sampled')}, purity(worst case) {r.get('purity_dual')}")


# --------------------------------------------------------------------------- C5
def task_c5(args):
    """Dependence on the size of the ambiguity set."""
    rows = []
    for sample in args.samples:
        af = args.aggregation_factor
        tag = f'c5_{sample}_af{af}'
        res = _solve(args.out, tag, sample=sample, aggregation_factor=af)
        params = make_params(sample=sample, aggregation_factor=af)
        data = rf.read_data(params)
        inst = baue_instanz(data[0], data[1], data[2], data[3], params)
        rows.append({
            'sample': sample,
            'delta_N': inst.zeit_diskret,
            'mu': list(inst.ret_time), 'mu_minus': list(inst.ret_time_minus),
            'mu_plus': list(inst.ret_time_plus),
            'sigma2_plus': list(inst.sigma2_plus),
            'schwank_var_global': inst.schwank_var_global,
            'result': _slim(res),
        })
    payload = {'task': 'C5', 'rows': rows}
    _dump(args.out, 'c5_ambiguity.json', payload)
    print(f"\n{'sample':>16} {'delta_N':>9} {'x^-':>9} {'x^+':>9} {'length':>8} "
          f"{'mass':>8} {'purity_wc':>10} {'purity_s':>9}")
    for r in rows:
        res = r['result']
        if res.get('interval_length') is None:
            print(f"{r['sample']:>16} {r['delta_N']:>9.6g} {'infeasible':>9}")
            continue
        print(f"{r['sample']:>16} {r['delta_N']:>9.6g} {res['x_lower']:>9.4f} {res['x_upper']:>9.4f} "
              f"{res['interval_length']:>8.4f} {res['collected_mass']:>8.4f} "
              f"{res['purity_dual'] or float('nan'):>10.5f} {res['purity_sampled']:>9.5f}")


# -------------------------------------------------------------------- C5 sweep
def task_c5_sweep(args):
    """Sweep over eps_ACN on regenerated data and locate the infeasibility threshold.

    The data for every eps_ACN is written by generate_data.py; the values inside
    [0.004, 0.0044] are interpolations between shipped data sets, all others are
    extrapolations of the calibrated surrogate of the map r_ACN -> mu_s(r_ACN),
    see the documentation of generate_data.py.
    """
    import generate_data as gen

    cal = gen.read_calibration()
    af = args.aggregation_factor
    rows = []

    def solve_eps(eps, name):
        gen.write_instance(cal, eps, name=name)
        res = _solve(args.out, f'c5_eps{eps:g}', sample=name, aggregation_factor=af,
                     time_limit=args.time_limit)
        params = make_params(sample=name, aggregation_factor=af)
        data = rf.read_data(params)
        inst = baue_instanz(data[0], data[1], data[2], data[3], params)
        # The model is always feasible: a window of length zero satisfies the purity
        # constraint as 0 >= 0. The quantity that decides whether the safe
        # approximation still says anything is the worst-case collected mass of the
        # desired species; if it is zero, the guaranteed purity is 0/0 and the
        # solution is vacuous.
        mass_wc = None
        if res.get('optimal_frac') is not None:
            lps = inner_lp.build_all(inst, 'lower', mip_indicator=True)
            i = params.wunschgroesse
            mass_wc = lps[i].solve(res['index_lower'], res['index_upper']) / inst.a_coeff(i)
        row = {'eps_ACN': eps, 'sample': name, 'delta_N': inst.zeit_diskret,
               'mu_minus': list(inst.ret_time_minus), 'mu_plus': list(inst.ret_time_plus),
               'mu_plus_minus_width': [p - m for p, m in zip(inst.ret_time_plus,
                                                             inst.ret_time_minus)],
               'sigma2_plus': list(inst.sigma2_plus),
               'feasible': res.get('optimal_frac') is not None,
               'collected_mass_worst_case': mass_wc,
               'informative': mass_wc is not None and mass_wc > 1e-6,
               'proven_optimal': res['stats'].get('proven_optimal'),
               'status': res['stats'].get('status_name'),
               'result': _slim(res)}
        rows.append(row)
        return row

    for eps in args.eps:
        r = solve_eps(eps, f'Gen_{eps:g}_E')
        res = r['result']
        print(f"  eps_ACN = {eps:g}: length {res.get('interval_length')}, "
              f"nominal mass {res.get('collected_mass')}, worst-case mass "
              f"{r['collected_mass_worst_case']}, purity(worst case) {res.get('purity_dual')}, "
              f"{'informative' if r['informative'] else 'VACUOUS (worst case collects nothing)'}")

    threshold = None
    if args.threshold:
        lo, hi = args.threshold_lower, args.threshold_upper
        while hi - lo > args.threshold_tol:
            mid = 0.5 * (lo + hi)
            r = solve_eps(mid, 'Bisect_E')
            print(f"  bisection eps_ACN = {mid:.6f}: length "
                  f"{r['result'].get('interval_length')}, worst-case mass "
                  f"{r['collected_mass_worst_case']}, "
                  f"{'informative' if r['informative'] else 'vacuous'}")
            if r['informative']:
                lo = mid
            else:
                hi = mid
        threshold = {'informative_up_to': lo, 'vacuous_from': hi,
                     'tolerance': args.threshold_tol}
        print(f"  the safe approximation still guarantees a positive worst-case yield for "
              f"eps_ACN <= {lo:.6f} and no longer for eps_ACN >= {hi:.6f}")
    payload = {'task': 'C5 sweep', 'aggregation_factor': af, 'rows': rows,
               'threshold': threshold, 'time_limit': args.time_limit}
    _dump(args.out, 'c5_eps_sweep.json', payload)


# ------------------------------------------------------------------------- eval
def task_eval(args):
    """Evaluate given fractionation intervals: collected mass, sampled purity and
    the certified enclosure of the worst-case purity over the ambiguity set."""
    af = args.aggregation_factor
    params = make_params(aggregation_factor=af)
    data = rf.read_data(params)
    inst = baue_instanz(data[0], data[1], data[2], data[3], params)
    t = np.asarray(inst.time_points)
    low = inner_lp.build_all(inst, 'lower')
    up = inner_lp.build_all(inst, 'upper')
    mip = inner_lp.build_all(inst, 'lower', mip_indicator=True)

    rows = []
    for spec in args.intervals:
        name, xl, xu = spec.split(':')
        xl, xu = float(xl), float(xu)
        lo = int(round((xl - t[0]) / inst.zeit_diskret))
        hi = int(round((xu - t[0]) / inst.zeit_diskret))
        frac = [1 if lo <= k <= hi else 0 for k in range(len(t))]
        vl = [l.solve(lo, hi) for l in low]
        vu = [u.solve(lo, hi) for u in up]
        vm = [l.solve(lo, hi) for l in mip]
        row = {
            'name': name,
            'x_lower_requested': xl, 'x_upper_requested': xu,
            'index_lower': lo, 'index_upper': hi,
            'x_lower': float(t[lo]), 'x_upper': float(t[hi]),
            'length': float(t[hi] - t[lo]),
            'collected_mass': aux.collected_mass(data, frac, params),
            'purity_sampled': aux.calculate_yieldpurity(data, frac, params),
            'val_lower': vl, 'val_upper': vu, 'val_mip_indicator': vm,
            'sum_val_lower': sum(vl), 'sum_val_upper': sum(vu), 'sum_val_mip_indicator': sum(vm),
            # worst-case collected mass of the desired species, i.e. min_P P([x^-,x^+]),
            # in contrast to 'collected_mass', which is the nominal one
            'collected_mass_worst_case_lower': vl[params.wunschgroesse] / inst.a_coeff(params.wunschgroesse),
            'collected_mass_worst_case_upper': vu[params.wunschgroesse] / inst.a_coeff(params.wunschgroesse),
            'collected_mass_worst_case_mip': vm[params.wunschgroesse] / inst.a_coeff(params.wunschgroesse),
            'purity_worst_case_lower': inner_lp.purity_from_values(inst, vl),
            'purity_worst_case_upper': inner_lp.purity_from_values(inst, vu),
            'purity_mip_indicator': inner_lp.purity_from_values(inst, vm),
        }
        rows.append(row)
        print(f"{name:>28}: [{row['x_lower']:.4f}, {row['x_upper']:.4f}] length {row['length']:.4f} "
              f"mass(nominal) {row['collected_mass']:.6f} "
              f"mass(worst case) in [{row['collected_mass_worst_case_lower']:.6f}, "
              f"{row['collected_mass_worst_case_upper']:.6f}] "
              f"purity(sampled) {row['purity_sampled']:.6f} "
              f"purity(worst case) in [{row['purity_worst_case_lower']:.6f}, "
              f"{row['purity_worst_case_upper']:.6f}] sum val in "
              f"[{row['sum_val_lower']:.6f}, {row['sum_val_upper']:.6f}]")
    _dump(args.out, args.name, {'task': 'eval', 'aggregation_factor': af,
                                'delta_N': inst.zeit_diskret, 'rows': rows})


# --------------------------------------------------------------------------- V1
def task_v1(args):
    """Moment bounds, variance bound and envelope maxima of every shipped data set."""
    rows = []
    for sample in args.samples:
        for af in args.factors:
            params = make_params(sample=sample, aggregation_factor=af)
            data = rf.read_data(params)
            inst = baue_instanz(data[0], data[1], data[2], data[3], params)
            ntp = 120000.
            rows.append({
                'sample': sample, 'aggregation_factor': af, 'delta_N': inst.zeit_diskret,
                'n_grid_points': len(inst.time_points),
                'T': [inst.time_points[0], inst.time_points[-1]],
                'q0_scaled': list(inst.q0),
                'mu': list(inst.ret_time), 'mu_minus': list(inst.ret_time_minus),
                'mu_plus': list(inst.ret_time_plus),
                'var_nom': list(inst.dict_var_nom), 'var_min': list(inst.dict_var_min),
                'var_max': list(inst.dict_var_max),
                'schwank_var_global': inst.schwank_var_global,
                'sigma2_plus': list(inst.sigma2_plus),
                'mu_squared_over_ntp': [m ** 2 / ntp for m in inst.ret_time],
                'sigma2_plus_over_mu_squared_over_ntp':
                    [s / (m ** 2 / ntp) for s, m in zip(inst.sigma2_plus, inst.ret_time)],
                'varianz_schranke': list(inst.varianz_schranke),
                'rho_max': [inst.rho_max(i) for i in inst.groessen],
                'envelope_mass': [aux.flaeche(inst.time_points, inst.schlauch_rtd[i])
                                  for i in inst.groessen],
                'a_coeff': [inst.a_coeff(i) for i in inst.groessen],
            })
            r = rows[-1]
            print(f"{sample:>8} af={af:<4} delta={r['delta_N']:.6g} n={r['n_grid_points']:<6} "
                  f"g={r['schwank_var_global']:.6f}")
            print(f"    mu_-  {['%.4f' % v for v in r['mu_minus']]}")
            print(f"    mu_+  {['%.4f' % v for v in r['mu_plus']]}")
            print(f"    s2+   {['%.5g' % v for v in r['sigma2_plus']]}")
            print(f"    s2+/(mu^2/NTP) {['%.6f' % v for v in r['sigma2_plus_over_mu_squared_over_ntp']]}")
            print(f"    rho_max {['%.4f' % v for v in r['rho_max']]}")
    _dump(args.out, 'v1_instance_data.json', {'task': 'V1', 'rows': rows})


# --------------------------------------------------------------------------- V2
def task_v2(args):
    """Rectangle rule versus exact envelope masses, and its effect on the model."""
    import envelope as env_mod

    af = args.aggregation_factor
    params = make_params(aggregation_factor=af)
    data = rf.read_data(params)
    inst = baue_instanz(data[0], data[1], data[2], data[3], params)
    rows = env_mod.comparison(inst)
    for r in rows:
        print(f"  species index {r['species_index']}: analytic fit residual "
              f"{r['fit_residual']:.2e}; sum(exact) {r['sum_exact']:.6f} vs sum(rect) "
              f"{r['sum_rectangle']:.6f}; cells where the rectangle rule is too small "
              f"{r['n_cells_rectangle_too_small']} (total deficit "
              f"{r['sum_positive_diff']:.3e}), too large {r['n_cells_rectangle_too_large']} "
              f"(total excess {r['sum_negative_diff']:.3e}); max relative deviation "
              f"{r['max_relative_diff']:.2e}")

    # effect on the adversarial values at the reported interval
    params_ex = make_params(aggregation_factor=af, exact_envelope_mass=True)
    inst_ex = baue_instanz(data[0], data[1], data[2], data[3], params_ex)
    t = np.asarray(inst.time_points)
    lo = int(round((args.x_lower - t[0]) / inst.zeit_diskret))
    hi = int(round((args.x_upper - t[0]) / inst.zeit_diskret))
    vals = {}
    for name, ins in (('rectangle', inst), ('exact', inst_ex)):
        lps = inner_lp.build_all(ins, 'lower', mip_indicator=True)
        vals[name] = [lp.solve(lo, hi) for lp in lps]
        print(f"  {name:>9}: val^s = {['%.6f' % v for v in vals[name]]}, "
              f"sum = {sum(vals[name]):.6f}, worst-case purity "
              f"{inner_lp.purity_from_values(ins, vals[name])}")

    payload = {'task': 'V2', 'aggregation_factor': af, 'delta_N': inst.zeit_diskret,
               'comparison': rows, 'interval': [float(t[lo]), float(t[hi])],
               'val_rectangle': vals['rectangle'], 'val_exact': vals['exact'],
               'sum_val_rectangle': sum(vals['rectangle']), 'sum_val_exact': sum(vals['exact']),
               'purity_rectangle': inner_lp.purity_from_values(inst, vals['rectangle']),
               'purity_exact': inner_lp.purity_from_values(inst_ex, vals['exact'])}

    if args.resolve:
        res = _solve(args.out, f'v2_exact_af{af}', aggregation_factor=af,
                     exact_envelope_mass=True, single_objective=True)
        payload['model_with_exact_masses'] = _slim(res)
        print(f"  model with exact envelope masses: interval "
              f"({res.get('x_lower')}, {res.get('x_upper')}), length {res.get('interval_length')}")
    _dump(args.out, 'v2_envelope_masses.json', payload)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', default='output', help="Output directory. Default 'output'.")
    sub = parser.add_subparsers(dest='task', required=True)

    sub.add_parser('c1').set_defaults(func=task_c1)

    p2 = sub.add_parser('c2')
    p2.add_argument('--time_limit', type=float, default=900.,
                    help="Time limit in seconds for the relaxed (enclosure) solves. Default 900.")
    p2.set_defaults(func=task_c2)

    p4 = sub.add_parser('c4')
    p4.add_argument('--aggregation_factor', type=int, default=1)
    p4.set_defaults(func=task_c4)

    p5 = sub.add_parser('c5')
    p5.add_argument('--aggregation_factor', type=int, default=1)
    p5.add_argument('--samples', nargs='+', default=['small', 'medium', 'large'])
    p5.set_defaults(func=task_c5)

    p5s = sub.add_parser('c5sweep')
    p5s.add_argument('--aggregation_factor', type=int, default=1)
    p5s.add_argument('--eps', type=float, nargs='+', default=[0.0021, 0.0042, 0.0063])
    p5s.add_argument('--threshold', action='store_true',
                     help="Locate by bisection the smallest eps_ACN at which the safe "
                          "approximation no longer guarantees a positive worst-case yield.")
    p5s.add_argument('--threshold_lower', type=float, default=0.0042,
                     help="Value of eps_ACN for which the safe approximation is informative.")
    p5s.add_argument('--threshold_upper', type=float, default=0.0063,
                     help="Value of eps_ACN for which it is vacuous.")
    p5s.add_argument('--threshold_tol', type=float, default=1e-4)
    p5s.add_argument('--time_limit', type=float, default=600.)
    p5s.set_defaults(func=task_c5_sweep)

    pe = sub.add_parser('eval')
    pe.add_argument('--aggregation_factor', type=int, default=1)
    pe.add_argument('--name', default='eval_intervals.json')
    pe.add_argument('--intervals', nargs='+',
                    default=['safe_approximation:3.2009:3.3813', 'nominal:3.1198:3.4760'],
                    help="Entries of the form name:x_lower:x_upper.")
    pe.set_defaults(func=task_eval)

    p2b = sub.add_parser('v2')
    p2b.add_argument('--aggregation_factor', type=int, default=1)
    p2b.add_argument('--x_lower', type=float, default=3.2009)
    p2b.add_argument('--x_upper', type=float, default=3.3813)
    p2b.add_argument('--resolve', action='store_true',
                     help="Also re-solve the MIP with the exact envelope masses.")
    p2b.set_defaults(func=task_v2)

    p1 = sub.add_parser('v1')
    p1.add_argument('--samples', nargs='+', default=['small', 'medium', 'large', 'long'])
    p1.add_argument('--factors', nargs='+', type=int, default=[1])
    p1.set_defaults(func=task_v1)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
