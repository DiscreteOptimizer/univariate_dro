"""Error bound of Theorem "inner_convergence" (equation "error_bound").

For every species ``s`` the bound

    Delta_N^s = |a^s| * ( (5 rho_max^s + 1) delta_N
                + 2 max{ 2 delta_N / ((mu_+^s - mu^s) + 2 delta_N),
                         ((2 Tbar + 1) delta_N + 5/4 delta_N^2)
                         / ((mu^s)^2 - mu_+^s mu_-^s + (2 Tbar + 1) delta_N + 5/4 delta_N^2) } )

with mu^s = (mu_-^s + mu_+^s)/2 and Tbar_s equal to the length of the
numerically relevant envelope support of species s bounds the difference
between the optimal value of the discretised inner problem and the optimal value
of the true inner problem. Summed over ``s`` it gives the amount by which the
right-hand side of the purity constraint (32b) has to be decreased to obtain the
upper bound val_N^up of Corollary "enclosure".

The coefficient a^s is the one that appears in the model, i.e.
``a_s[i] * q0[i]`` including the mass scaling ``q_factor``; see
``instanz.Instanz.a_coeff``.
"""

from typing import Dict, List

from instanz import Instanz


def delta_bound(inst: Instanz, i: int) -> Dict[str, float]:
    """Delta_N^s of species ``i`` together with its ingredients."""
    delta = inst.zeit_diskret
    t_bar = inst.t_bar(i)
    rho_max = inst.rho_max(i)
    mu_minus = inst.ret_time_minus[i]
    mu_plus = inst.ret_time_plus[i]
    mu = 0.5 * (mu_minus + mu_plus)
    a = abs(inst.a_coeff(i))

    term_env = (5. * rho_max + 1.) * delta

    # first moment term
    term_mean = 2. * delta / ((mu_plus - mu) + 2. * delta)

    # second moment term
    num = (2. * t_bar + 1.) * delta + 1.25 * delta ** 2
    den = (mu ** 2 - mu_plus * mu_minus) + num
    term_var = num / den

    bound = a * (term_env + 2. * max(term_mean, term_var))
    return {
        'species_index': i,
        'delta_N': delta,
        'T_bar': t_bar,
        'rho_max': rho_max,
        'mu_minus': mu_minus,
        'mu_plus': mu_plus,
        'mu': mu,
        'abs_a': a,
        'term_envelope': term_env,
        'term_first_moment': term_mean,
        'term_second_moment': term_var,
        'eta_N_admissible': delta ** 2 / (2. * rho_max * (t_bar ** 3 + 1.)),
        'Delta_N': bound,
    }


def delta_bounds(inst: Instanz) -> List[Dict[str, float]]:
    """Delta_N^s for all species."""
    return [delta_bound(inst, i) for i in inst.groessen]


def total_delta(inst: Instanz) -> float:
    """sum_s Delta_N^s, i.e. the offset of the right-hand side of (32b)."""
    return sum(d['Delta_N'] for d in delta_bounds(inst))


def _main():
    import argparse
    import json

    import run_funktionen as rf
    from instanz import baue_instanz
    from params import Params

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--aggregation_factor", type=int, default=1)
    parser.add_argument("--sample", type=str, default='long')
    parser.add_argument("--reinheit", type=float, default=.95)
    parser.add_argument("--wunschgroesse", type=int, default=2)
    parser.add_argument("--json", type=str, default=None)
    args = parser.parse_args()

    params = Params(aggregation_factor=args.aggregation_factor, reinheit=args.reinheit,
                    wunschgroesse=args.wunschgroesse, sample=args.sample, fix=False,
                    fix_lower=0, fix_upper=0, nominal=False, plot=False)
    data = rf.read_data(params)
    inst = baue_instanz(data[0], data[1], data[2], data[3], params)
    bounds = delta_bounds(inst)
    for b in bounds:
        print(f"s index {b['species_index']}: |a^s|={b['abs_a']:.6f} rho_max={b['rho_max']:.4f} "
              f"env={b['term_envelope']:.6e} mean={b['term_first_moment']:.6e} "
              f"var={b['term_second_moment']:.6e} Delta_N={b['Delta_N']:.6e}")
    print(f"delta_N = {inst.zeit_diskret:.6g}")
    print(f"sum_s Delta_N^s = {sum(b['Delta_N'] for b in bounds):.6f}")
    if args.json:
        with open(args.json, 'w') as f:
            json.dump({'delta_N': inst.zeit_diskret,
                       'T_bar': [inst.t_bar(i) for i in inst.groessen],
                       'bounds': bounds,
                       'total': sum(b['Delta_N'] for b in bounds)}, f, indent=2)


if __name__ == "__main__":
    _main()
