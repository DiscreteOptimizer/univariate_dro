# general imports
import time
import argparse

# auxiliary functions
import hilfsfunktionen as aux
from params import Params

# sub modules
import einlesen as reader
import model_2 as m2
import my_plot as plotter


def read_data(params: Params):
    """Read, aggregate and trim the input data of one instance."""
    # Einlesen
    daten, name = reader.rd4(params)

    # hacky nominal optimization cheat
    if params.nominal:
        daten[1] = daten[2] = daten[0]

    # print data paths
    print(f'data_nom: {daten[0]}')
    print(f'data_min: {daten[1]}')
    print(f'data_max: {daten[2]}')

    """
        The input consists of three distributions per particle,
        nominal (expected traverse),
        minimal (fast traverse),
        maximal (slow traverse).
    """

    start = time.time()

    # Load input data
    data = list(reader.inputdata(daten[0], daten[1], daten[2]))
    print(f"Time to read input: {time.time() - start:.2f}s")

    time_points = aux.aggregate_matrix_index(data[0], params.aggregation_factor)
    matrices = [aux.aggregate_matrix(data[i], params.aggregation_factor, params) for i in (1, 2, 3)]
    data[0] = time_points
    data[1], data[2], data[3] = matrices

    data = aux.remove_zeros(data)
    return data


def run(params: Params):
    print(f'\n---Start---')

    data = read_data(params)

    # Solve model
    # opt contains the fract vector
    opt = m2.solve_dro_model(data[0], data[1], data[2], data[3], params)

    if opt['optimal_frac'] is None:
        print("[INFEASIBLE] No fractionation interval exists for this instance.")
        return opt

    # calculate purity
    purity = aux.calculate_yieldpurity(data, opt['optimal_frac'], params)
    opt['purity_sampled'] = purity
    # collected mass of the desired species, relative to its total mass
    opt['collected_mass'] = aux.collected_mass(data, opt['optimal_frac'], params)

    # some cln output
    if params.reinheit > purity:
        print(f"[FAILED] Best purity reached: Purity {purity:.4f} < Target {params.reinheit:.4f}.")
    else:
        print(f"[SUCCESS] Target (nominal) purity reached: Purity {purity:.4f} >= Target {params.reinheit:.4f}.")
    print(f"Collected mass (relative to total mass of desired species): {opt['collected_mass']:.6f}")

    if params.stats_file:
        # re-write, now including the sampled purity and the collected mass
        m2.write_stats(opt, params.stats_file)

    # generate plot
    if params.plot:
        plotter.plot(data, opt["optimal_frac"], params)

    return opt

# ----------------------------
# main function with argparse
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Robust Chromatography with DRO."
    )

    parser.add_argument("--aggregation_factor", type=int, default=1,
                        help="Aggregation parameter. Default 1.")
    parser.add_argument("--reinheit", type=float, default=.95,
                        help="Desired purity. Default 0.95.")
    parser.add_argument("--wunschgroesse", type=int, default=2,
                        help="Desired particle size. Default 2.")
    parser.add_argument("--fix_lower", type=int, default=299,
                        help="Fix lower fract bound. Default 299.")
    parser.add_argument("--fix_upper", type=int, default=653,
                        help="Fix upper fract bound. Default 653.")
    parser.add_argument("--sample", type=str, default='long',
                        help="Sample. 'small', 'medium', 'large', 'long' (default), or the "
                             "name of a data directory below 'daten' without the _min/_max suffix.")
    parser.add_argument("-f", "--fix", action="store_true",
                        help="Fix solution to fix_lower, fix_upper. Check purity for this solution.")
    parser.add_argument("-n", "--nominal", action="store_true",
                        help="Solve nominal program instead. Hacky.")
    # --- added in 0.2.0 ---
    parser.add_argument("--purity_rhs", type=float, default=0.0,
                        help="Right-hand side of the purity constraint (32b). 0.0 (default) is "
                             "the safe approximation; pass -sum_s Delta_N^s (see error_bound.py) "
                             "to obtain the upper bound of the certified enclosure.")
    parser.add_argument("--exact_envelope_mass", action="store_true",
                        help="Use the exact integral of the envelope over every grid cell "
                             "instead of the rectangle rule (aggregation_factor 1 only).")
    parser.add_argument("--no_second_moment", action="store_true",
                        help="Drop the relaxed second moment constraint by fixing its dual "
                             "variable to zero.")
    parser.add_argument("--log_file", type=str, default=None,
                        help="Write the Gurobi log to this file.")
    parser.add_argument("--stats_file", type=str, default=None,
                        help="Write model size, solver statistics and solution as JSON to this file.")
    parser.add_argument("--no_plot", action="store_true",
                        help="Do not produce the matplotlib figure.")

    args = parser.parse_args()

    params = Params(
        aggregation_factor=args.aggregation_factor,
        reinheit=args.reinheit,
        wunschgroesse=args.wunschgroesse,
        sample=args.sample,
        fix=args.fix,
        fix_lower=args.fix_lower,
        fix_upper=args.fix_upper,
        nominal=args.nominal,
        purity_rhs=args.purity_rhs,
        no_second_moment=args.no_second_moment,
        exact_envelope_mass=args.exact_envelope_mass,
        log_file=args.log_file,
        stats_file=args.stats_file,
        plot=not args.no_plot,
    )

    run(params)

if __name__ == "__main__":
    main()
