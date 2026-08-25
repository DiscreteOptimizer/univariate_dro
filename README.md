# univariate_dro

Safe Mixed-Integer Approximation for Non-Convex Distributional Robustness with Univariate Indicator Functions

## Short Description

This repository contains the code to reproduce the computational results of the paper

> Dienstbier, J.; Liers, F.; Rolfes, J.; Rösel, F. (2026):
> *Safe Mixed-Integer Approximation for Non-Convex Distributional Robustness with Univariate Indicator Functions*,
> submitted to the *Journal of Global Optimization*.

---

## Requirements

- **Python:** 3.12.
- **Operating system:** tested under Ubuntu 24.04.
- For the code to run, Gurobi must be installed and a valid licence must exist.
- Further dependencies: see 'requirements.txt'.

---

## Installation

The recommended setup uses a virtual environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The directory daten and all provided data files must be located in the working directory from which the code is executed.

---

## Execution

```bash
# Robust optimization
python3.12 run_funktionen.py

# Nominal optimization
python3.12 run_funktionen.py -n

# Evaluation of fixed fractionation times
python3.12 run_funktionen.py -f --fix_lower <number1> --fix_upper <number2>

# Show all available options
python3.12 run_funktionen.py -h
```

The most important further options of `run_funktionen.py` are

| option | effect |
|---|---|
| `--aggregation_factor <k>` | coarsen the time grid by a factor `k`, i.e. `delta_N = k * 1e-4` for the default data set |
| `--purity_rhs <r>` | right-hand side of the purity constraint (32b). `0.0` is the safe approximation; `-sum_s Delta_N^s` gives the upper bound of the certified enclosure, see `error_bound.py` |
| `--no_second_moment` | drop the relaxed second-moment constraint by fixing its dual variable to zero |
| `--exact_envelope_mass` | integrate the envelope over each grid cell exactly instead of using the rectangle rule, see `envelope.py` |
| `--refinement_factor <r>` | refine the time grid by a zero-order hold, i.e. `delta_N = 1e-4 / r` for the default data set. Requires `--aggregation_factor 1` |
| `--refine_recompute` | with the above: recompute the moment bounds and the envelope on the refined grid instead of taking them from the data grid. Changes the ambiguity set slightly; for diagnostics only |
| `--stats_file <path>` | write model size, solver statistics and the solution as JSON |
| `--log_file <path>` | write the Gurobi log to a file |
| `--no_plot` | do not produce the matplotlib figure |

---

## Modules

| file | content |
|---|---|
| `run_funktionen.py` | entry point: reads the data, solves the model, evaluates and plots the solution |
| `einlesen.py` | reading of the data directories |
| `hilfsfunktionen.py` | grid aggregation, envelope, areas, sampled purity, collected mass |
| `instanz.py` | pre-processing of one instance: particle masses, normalised densities, moment bounds, variance bound, envelope |
| `model_2.py` | the mixed-integer program and its solution statistics |
| `my_plot.py` | chromatogram figure of a solution |
| `params.py` | parameters of one run |
| `error_bound.py` | the a priori error bound `Delta_N^s` of the convergence theorem |
| `inner_lp.py` | linear programs that bracket the adversarial problem of one species from below and from above |
| `reference.py` | reference solution of the original problem by enumeration of the two outer variables |
| `worst_case.py` | worst-case measures of the safe approximation and the corresponding figure |
| `envelope.py` | exact envelope masses per grid cell (comparison with the rectangle rule) |
| `generate_data.py` | regeneration of the residence time distributions for other values of `eps_ACN` |
| `experiments.py` | driver that produces all reported numbers and writes them as JSON |
| `check_refinement.py` | consistency checks of the grid refinement |
| `slide_plot.py` | chromatogram figure with the worst-case densities, one vector PDF per fractionation window |

---

## Reproducing the reported numbers

```bash
# model sizes, running times and B&B nodes for all discretisations
python3.12 experiments.py --out results c1

# certified enclosure of the optimal value
python3.12 experiments.py --out results c2 --time_limit 900

# effect of the second-moment constraint
python3.12 experiments.py --out results c4

# dependence on the size of the ambiguity set
python3.12 generate_data.py --validate
python3.12 generate_data.py --eps 0.0021 0.0063
python3.12 experiments.py --out results c5 --samples Gen_0.0021_E long Gen_0.0063_E

# moment and variance bounds of every data set, envelope masses
python3.12 experiments.py --out results v1
python3.12 experiments.py --out results v2 --resolve

# evaluation of given fractionation intervals
python3.12 experiments.py --out results eval \
    --intervals "safe_approximation:3.2009:3.3813" "nominal:3.1198:3.4760"

# finer discretisations: delta_N = 1e-4 / r
python3.12 check_refinement.py --refinement 2 5 10 --json results/r1_checks.json
python3.12 experiments.py --out results r1 --refinement 1 2 5 10 --time_limit 2400

# worst-case densities of least curvature, one panel per window
python3.12 slide_plot.py --window robust:3.2009:3.3813 --objective curvature \
    --out results/robust.pdf --json results/worst_case_robust.json

# reference solution by enumeration
python3.12 reference.py --json results/reference.json

# worst-case measures
python3.12 worst_case.py --out results/worst_case_measures.pdf --json results/worst_case.json
```

Every sub-command of `experiments.py` writes one JSON file with all numbers and one
Gurobi log per solve into `--out`, so that each number can be traced back to a raw
solver log.

---

## Input data

The input data is located in the directory "daten". It is expected that it contains three folders ("nom", "min", "max"), each containing %num-particles many .txt files with lines in the format
time_point<blank>particle_density

The shipped data sets are

| directory | `eps_ACN` | grid width |
|---|---|---|
| `Small_E` | 0.0040 | 1e-3 |
| `Medium_E` | 0.0042 | 1e-3 |
| `Large_E` | 0.0044 | 1e-3 |
| `Long_E` (default) | 0.0042 | 1e-4 |

`Long_E` contains the same distributions as `Medium_E` on a ten times finer grid.
The parameter `sample` accepts the values "small", "medium", "large", "long"
(default), and additionally the name of any other data directory below `daten`,
which is what `generate_data.py` writes.

---

## Output data

If the program terminates successfully, it produces command-line output and a matplotlib figure plot.pdf, which is written to the (automatically generated) directory ./output.

---

## License

The license under which the code is provided is specified in the file LICENSE.

---

## Citation

If you want to cite this code, please cite the above-mentioned article.
