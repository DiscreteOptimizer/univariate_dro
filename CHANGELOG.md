# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [0.1.0] - 2025-12-19

### Added
- Initial public release of the code.
- Implementation of the mixed-integer approximation approach for non-convex distributionally robust fractionation.
- Scripts to reproduce the computational results reported in the associated article.
- Documentation and example usage via `README.md`.

### Notes
- This release reproduces the computational results of the manuscript
  *“Safe Mixed-Integer Approximation for Non-Convex Distributional Robustness with Univariate Indicator Functions”*,
  submitted to the *European Journal of Operational Research*.

## [0.1.1] - 2026-01-18

### Minor Fixes
- Added Long instance set (default), with a minimum time discretization of .0001 instead of .001.
- Set default (integer and constraint) feasibility tolerances to 1e-09, default MIPGap to 0., and NumericFocus to 3.

## [0.2.0] - 2026-08-17

This release adds the computations requested in the first review round of the
manuscript (now under review at the *Journal of Global Optimization*). The model
itself is unchanged: every new behaviour is behind a command line switch, and the
default run reproduces the results of version 0.1.1 exactly.

### Added
- `instanz.py`: the pre-processing of an instance (particle masses, normalised
  densities, moment bounds, variance bound, envelope) was moved out of
  `model_2.solve_dro_model` into `baue_instanz`, so that the analysis tools work
  on exactly the same data as the model. The arithmetic is unchanged.
- `error_bound.py`: the a priori error bound `Delta_N^s` of the convergence
  theorem.
- `inner_lp.py`: two linear programs per species that bracket the value of the
  adversarial problem from below (the re-dualised discretised problem, which is
  the relaxation the safe approximation is built on) and from above (a
  restriction to measures with cell-wise constant density). The lower bound LP
  can use either the indicator of the convergence proof or the weaker indicator
  of the mixed-integer model.
- `reference.py`: reference solution of the original problem by enumeration of
  the two outer variables, with a certified enclosure of its optimal value.
- `worst_case.py`: the worst-case measures of the safe approximation and a figure
  showing that their density stays below the envelope.
- `envelope.py`: exact envelope mass of every grid cell, for comparison with the
  rectangle rule used in the model.
- `generate_data.py`: regeneration of the residence time distributions for other
  values of `eps_ACN`, calibrated on the shipped data sets, including a
  validation against them.
- `experiments.py`: driver that produces all reported numbers and writes them as
  JSON together with one Gurobi log per solve.
- New options of `run_funktionen.py`: `--purity_rhs`, `--no_second_moment`,
  `--exact_envelope_mass`, `--stats_file`, `--log_file`, `--no_plot`; `--sample`
  now also accepts the name of an arbitrary data directory below `daten`.
- New options of `Params`: `purity_rhs`, `no_second_moment`,
  `exact_envelope_mass`, `single_objective`, `time_limit`, `log_file`,
  `stats_file`, `plot`, `keep_model`. All default to the behaviour of 0.1.1.
- `model_2.solve_dro_model` now returns the model size, the solver statistics
  (running time, B&B nodes, iterations, dual bound, status), the fractionation
  interval in grid indices and in minutes, the per-species values of the purity
  constraint and the worst-case purity.
- `hilfsfunktionen.collected_mass`: collected mass of the desired species inside
  the fractionation window, relative to its total mass.

### Fixed
- The fractionation vector is now rounded instead of truncated when it is read
  from the solution (`int(round(x))` instead of `int(x)`), so that an integer
  value returned as `1 - 1e-11` is no longer read as `0`.
- `solve_dro_model` no longer raises when a model has no solution; it reports the
  status and returns a result without a fractionation interval.
- Multi-objective models do not expose `MIPGap`; querying solver attributes is
  now guarded.

### Changed
- The article is now under review at the *Journal of Global Optimization*;
  `README.md` and `CITATION.cff` were updated accordingly.
