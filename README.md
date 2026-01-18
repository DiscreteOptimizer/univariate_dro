# univariate_dro

Safe Mixed-Integer Approximation for Non-Convex Distributional Robustness with Univariate Indicator Functions

## Short Description

This repository contains the code to reproduce the computational results of the paper

> Dienstbier, J.; Liers, F.; Rolfes, J.; Rösel, F. (2025):
> *Safe Mixed-Integer Approximation for Non-Convex Distributional Robustness with Univariate Indicator Functions*,
> submitted to the *European Journal of Operational Research*.

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

---

## Input data

The input data is located in the directory "daten". It is expected that it contains three folders ("nom", "min", "max"), each containing %num-particles many .txt files with lines in the format
time_point<blank>particle_density

Which folders are used is currently hard-coded and can be steered with the input parameter "sample", which accepts the values "small", "medium", "large", and "long" (default). If you want to solve other instances, you have to adapt the code.

---

## Output data

If the program terminates successfully, it produces command-line output and a matplotlib figure plot.pdf, which is written to the (automatically generated) directory ./output.

---

## License

The license under which the code is provided is specified in the file LICENSE.

---

## Citation

If you want to cite this code, please cite the above-mentioned article.

