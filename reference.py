"""Reference solution of Problem (35) ("Prob: DRO_chromatography") by enumeration.

The outer problem has only the two variables ``x^-`` and ``x^+``, so it can be
solved by enumeration without using the safe approximation: for every candidate
pair the four adversarial problems are solved as linear programs, see
``inner_lp.py``.

Two tests are used for every candidate pair, both restricted to grid points:

certified
    ``sum_s LowerLP^s >= 0``. Since ``LowerLP^s <= val^s``, a pair that passes is
    feasible for the true Problem (35). The best length over all passing pairs is
    therefore a lower bound on the optimal value ``val^*``.

optimistic
    ``sum_s UpperLP^s >= 0`` with ``cover=True``. Since ``UpperLP^s >= val^s``, a
    pair that fails is infeasible for (35), and with ``cover=True`` this holds for
    every real interval within one grid cell of the tested pair. The best length
    over all passing pairs is therefore an upper bound on ``val^*``.

The optimal value of (35) lies between the two, which is the enclosure reported
in the manuscript's reference table.

The scan is organised in two phases. Phase 1 covers the whole domain on a coarse
evaluation grid; because the optimistic test is valid on any grid, this phase
already excludes most of the domain rigorously. Phase 2 repeats both tests on the
fine grid, restricted to the neighbourhood of the pairs that survived phase 1.
"""

import argparse
import json
import time

import numpy as np

import inner_lp
import run_funktionen as rf
from instanz import baue_instanz
from params import Params


class Scanner:
    """Evaluates ``sum_s val^s`` for intervals given by grid indices."""

    def __init__(self, inst, kind, cover=False):
        self.inst = inst
        self.kind = kind
        self.lps = inner_lp.build_all(inst, kind)
        if cover:
            for lp in self.lps:
                lp.cover = True
        # contaminants first: they contribute negatively, so the partial sum plus
        # the largest possible remaining contribution is an upper bound on the
        # total and allows an early exit
        order = sorted(inst.groessen, key=lambda i: inst.a_coeff(i))
        self.order = order
        self.rest_ub = []
        for pos in range(len(order)):
            self.rest_ub.append(sum(max(inst.a_coeff(i), 0.) for i in order[pos + 1:]))
        self.n_lp = 0

    def evaluate(self, lo, hi, prune=True):
        """Return ``(feasible, total, values)``; ``values`` may be partial when
        the evaluation was cut short by the early exit."""
        values = {}
        total = 0.
        for pos, i in enumerate(self.order):
            v = self.lps[i].solve(lo, hi)
            self.n_lp += 1
            values[i] = v
            total += v
            if prune and total + self.rest_ub[pos] < 0.:
                return False, total + self.rest_ub[pos], values
        return total >= 0., total, values


def scan_region(scanner, lo_range, hi_range, min_len, verbose_every=0, t=None):
    """All pairs in ``lo_range x hi_range`` with ``hi - lo >= min_len`` that pass.

    Pairs are enumerated with decreasing length, and for every length the scan
    stops as soon as one pair passes: the first hit is an optimal pair on the
    candidate grid, so no shorter length has to be looked at.
    """
    best = None
    passing = []
    max_len = hi_range[-1] - lo_range[0]
    lengths = range(max_len, min_len - 1, -1)
    start = time.time()
    for k, length in enumerate(lengths):
        hits = []
        for lo in lo_range:
            hi = lo + length
            if hi not in hi_range:
                continue
            ok, total, values = scanner.evaluate(lo, hi)
            if ok:
                hits.append((lo, hi, total, values))
        if verbose_every and k % verbose_every == 0:
            print(f"    length {length} cells: {len(hits)} hits, {scanner.n_lp} LPs, "
                  f"{time.time() - start:.1f}s", flush=True)
        if hits:
            passing = hits
            best = hits[0]
            break
    return best, passing


def full_scan(scanner, n_cells, min_len, lo_lim=None, verbose_every=200):
    """Scan the whole domain (all pairs) and return every passing pair."""
    lo_lim = lo_lim if lo_lim is not None else (0, n_cells)
    passing = []
    start = time.time()
    for lo in range(lo_lim[0], lo_lim[1] + 1):
        for hi in range(lo + min_len, n_cells + 1):
            ok, total, values = scanner.evaluate(lo, hi)
            if ok:
                passing.append((lo, hi, total, values))
        if verbose_every and lo % verbose_every == 0:
            print(f"    lo={lo}: {len(passing)} passing pairs so far, {scanner.n_lp} LPs, "
                  f"{time.time() - start:.1f}s", flush=True)
    return passing


def _instance(sample, af, reinheit, wunschgroesse):
    params = Params(aggregation_factor=af, reinheit=reinheit, wunschgroesse=wunschgroesse,
                    sample=sample, fix=False, fix_lower=0, fix_upper=0, nominal=False,
                    plot=False)
    data = rf.read_data(params)
    inst = baue_instanz(data[0], data[1], data[2], data[3], params)
    return data, inst, params


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sample", type=str, default='long')
    parser.add_argument("--reinheit", type=float, default=.95)
    parser.add_argument("--wunschgroesse", type=int, default=2)
    parser.add_argument("--coarse_factor", type=int, default=10,
                        help="Aggregation factor of the coarse evaluation grid. Default 10, "
                             "i.e. delta = 0.001.")
    parser.add_argument("--fine_factor", type=int, default=1,
                        help="Aggregation factor of the fine evaluation grid. Default 1, "
                             "i.e. delta = 0.0001.")
    parser.add_argument("--coarse_min_len", type=float, default=0.05,
                        help="Shortest interval considered in the coarse full scan, in minutes.")
    parser.add_argument("--margin", type=float, default=0.002,
                        help="The fine scan covers the bounding box of the coarse passing pairs, "
                             "enlarged by this margin in minutes.")
    parser.add_argument("--json", type=str, default=None)
    args = parser.parse_args()

    out = {'settings': vars(args), 'phases': {}}

    # ------------------------------------------------------------------ phase 1
    print("=== phase 1: coarse full scan over the whole domain")
    data_c, inst_c, _ = _instance(args.sample, args.coarse_factor, args.reinheit, args.wunschgroesse)
    t_c = np.asarray(inst_c.time_points)
    d_c = inst_c.zeit_diskret
    min_len_c = int(round(args.coarse_min_len / d_c))
    print(f"    grid: {len(t_c)} points, delta = {d_c:.6g}, "
          f"T = [{t_c[0]:.5f}, {t_c[-1]:.5f}], min length {min_len_c} cells")

    opt_c = Scanner(inst_c, 'upper', cover=True)
    start = time.time()
    passing_c = full_scan(opt_c, inst_c.anzahl_prozess, min_len_c)
    wall_opt_c = time.time() - start
    print(f"    optimistic test: {len(passing_c)} passing pairs, {opt_c.n_lp} LPs, {wall_opt_c:.1f}s")
    if not passing_c:
        print("    no pair passes the optimistic test: Problem (35) is infeasible "
              "for intervals of the considered lengths.")
        out['phases']['coarse_optimistic'] = {'passing': 0, 'wall_time': wall_opt_c}
        _dump(out, args.json)
        return

    best_opt_c = max(passing_c, key=lambda r: r[1] - r[0])
    print(f"    longest optimistic pair: ({t_c[best_opt_c[0]]:.5f}, {t_c[best_opt_c[1]]:.5f}), "
          f"length {(best_opt_c[1] - best_opt_c[0]) * d_c:.5f}")
    out['phases']['coarse_optimistic'] = {
        'delta_N': d_c, 'passing': len(passing_c), 'wall_time': wall_opt_c, 'n_lp': opt_c.n_lp,
        'best': _pair(t_c, best_opt_c, d_c),
    }

    # The certified test only has to be applied to the pairs that survived the
    # optimistic one, since LowerLP <= val^s <= UpperLP for every pair. Testing
    # them by decreasing length and stopping at the first success gives the
    # longest certified pair on the coarse candidate grid.
    cert_c = Scanner(inst_c, 'lower')
    start = time.time()
    best_cert_c = None
    for lo, hi, _, _ in sorted(passing_c, key=lambda r: r[0] - r[1]):
        ok, total, values = cert_c.evaluate(lo, hi)
        if ok:
            best_cert_c = (lo, hi, total, values)
            break
    wall_cert_c = time.time() - start
    if best_cert_c is None:
        print(f"    certified test: no pair passes on the coarse grid ({wall_cert_c:.1f}s)")
        out['phases']['coarse_certified'] = {'delta_N': d_c, 'best': None, 'wall_time': wall_cert_c}
    else:
        print(f"    longest certified pair: ({t_c[best_cert_c[0]]:.5f}, {t_c[best_cert_c[1]]:.5f}), "
              f"length {(best_cert_c[1] - best_cert_c[0]) * d_c:.5f} ({wall_cert_c:.1f}s)")
        out['phases']['coarse_certified'] = {'delta_N': d_c, 'best': _pair(t_c, best_cert_c, d_c),
                                            'wall_time': wall_cert_c, 'n_lp': cert_c.n_lp}

    # region that still has to be searched on the fine grid: every pair that
    # passed the optimistic test and is at least as long as the certified one
    len_cert = (best_cert_c[1] - best_cert_c[0]) if best_cert_c else min_len_c
    relevant = [r for r in passing_c if r[1] - r[0] >= len_cert]
    lo_min = min(r[0] for r in relevant)
    lo_max = max(r[0] for r in relevant)
    hi_min = min(r[1] for r in relevant)
    hi_max = max(r[1] for r in relevant)
    print(f"    {len(relevant)} pairs are at least as long as the certified one; "
          f"their bounding box is x^- in [{t_c[lo_min]:.5f}, {t_c[lo_max]:.5f}], "
          f"x^+ in [{t_c[hi_min]:.5f}, {t_c[hi_max]:.5f}]")
    out['phases']['coarse_optimistic']['box'] = [t_c[lo_min], t_c[lo_max], t_c[hi_min], t_c[hi_max]]
    out['phases']['coarse_optimistic']['n_relevant'] = len(relevant)

    # ------------------------------------------------------------------ phase 2
    print("=== phase 2: fine scan in the neighbourhood of the coarse passing pairs")
    data_f, inst_f, _ = _instance(args.sample, args.fine_factor, args.reinheit, args.wunschgroesse)
    t_f = np.asarray(inst_f.time_points)
    d_f = inst_f.zeit_diskret
    print(f"    grid: {len(t_f)} points, delta = {d_f:.6g}")

    def fine_index(x):
        return int(round((x - t_f[0]) / d_f))

    margin = args.margin
    lo_lo = max(fine_index(t_c[lo_min] - margin), 0)
    lo_hi = min(fine_index(t_c[lo_max] + margin), inst_f.anzahl_prozess)
    hi_lo = max(fine_index(t_c[hi_min] - margin), 0)
    hi_hi = min(fine_index(t_c[hi_max] + margin), inst_f.anzahl_prozess)
    # the fine optimum cannot be shorter than the coarse certified one
    min_len_f = int(round(((best_cert_c[1] - best_cert_c[0]) * d_c) / d_f)) if best_cert_c else \
        int(round(args.coarse_min_len / d_f))
    print(f"    searched region: x^- in [{t_f[lo_lo]:.5f}, {t_f[lo_hi]:.5f}], "
          f"x^+ in [{t_f[hi_lo]:.5f}, {t_f[hi_hi]:.5f}], min length {min_len_f} cells")

    res_fine = {}
    for kind, cover in (('upper', True), ('lower', False)):
        sc = Scanner(inst_f, kind, cover=cover)
        start = time.time()
        best, passing = scan_region(sc, range(lo_lo, lo_hi + 1), range(hi_lo, hi_hi + 1),
                                    min_len_f, verbose_every=20, t=t_f)
        wall = time.time() - start
        name = 'fine_optimistic' if kind == 'upper' else 'fine_certified'
        if best is None:
            print(f"    {name}: no pair passes ({wall:.1f}s, {sc.n_lp} LPs)")
            res_fine[name] = {'delta_N': d_f, 'best': None, 'wall_time': wall, 'n_lp': sc.n_lp}
        else:
            lo, hi, total, values = best
            print(f"    {name}: ({t_f[lo]:.5f}, {t_f[hi]:.5f}), length {(hi - lo) * d_f:.5f}, "
                  f"sum val = {total:.6f} ({wall:.1f}s, {sc.n_lp} LPs, {len(passing)} pairs at "
                  f"this length)")
            res_fine[name] = {'delta_N': d_f, 'best': _pair(t_f, best, d_f), 'wall_time': wall,
                              'n_lp': sc.n_lp, 'n_optimal_pairs': len(passing),
                              'all_optimal_pairs': [_pair(t_f, r, d_f) for r in passing]}
    out['phases'].update(res_fine)

    # ---------------------------------------------------------------- summary
    lower = res_fine.get('fine_certified', {}).get('best')
    upper = res_fine.get('fine_optimistic', {}).get('best')
    if lower and upper:
        print(f"=== enclosure of the optimal value of (35): "
              f"[{lower['length']:.5f}, {upper['length']:.5f}] minutes")
        out['enclosure'] = {'lower': lower['length'], 'upper': upper['length'],
                            'x_lower': lower['x_lower'], 'x_upper': lower['x_upper']}
    _dump(out, args.json)


def _pair(t, rec, d):
    lo, hi, total, values = rec
    return {'index_lower': int(lo), 'index_upper': int(hi),
            'x_lower': float(t[lo]), 'x_upper': float(t[hi]),
            'length': float((hi - lo) * d), 'sum_val': float(total),
            'values': {str(k): float(v) for k, v in values.items()}}


def _dump(out, path):
    if path:
        with open(path, 'w') as f:
            json.dump(out, f, indent=2, default=float)
        print(f"written: {path}")


if __name__ == "__main__":
    main()
