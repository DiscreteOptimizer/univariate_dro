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

def run(params: Params):
    print(f'\n---Start---')
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
    # scaling ist am anfang auf 1.0 - passiert da einfach nichts?
    data[0] = time_points
    data[1], data[2], data[3] = matrices

    data = aux.remove_zeros(data)
    # Solve model
    # opt contains the fract vector
    opt = m2.solve_dro_model(data[0], data[1], data[2], data[3], params)
    # calculate purity
    purity = aux.calculate_yieldpurity(data, opt['optimal_frac'], params)

    # some cln output
    if params.reinheit > purity:
        print(f"[FAILED] Best purity reached: Purity {purity:.4f} < Target {params.reinheit:.4f}.")
    else:
        print(f"[SUCCESS] Target (nominal) purity reached: Purity {purity:.4f} >= Target {params.reinheit:.4f}.")

    # generate plot
    plotter.plot(data, opt["optimal_frac"], params)

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
                        help="Sample. Default 'medium'.")
    parser.add_argument("-f", "--fix", action="store_true",
                        help="Fix solution to fix_lower, fix_upper. Check purity for this solution.")
    parser.add_argument("-n", "--nominal", action="store_true",
                        help="Solve nominal program instead. Hacky.")

    args = parser.parse_args()

    params = Params(
        aggregation_factor=args.aggregation_factor,
        reinheit=args.reinheit,
        wunschgroesse=args.wunschgroesse,
        sample=args.sample,
        fix=args.fix,
        fix_lower=args.fix_lower,
        fix_upper=args.fix_upper,
        nominal=args.nominal
    )

    run(params)

if __name__ == "__main__":
    main()

