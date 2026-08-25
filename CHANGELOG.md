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

## [0.3.0] - 2026-08-25

Round 2 of the JOGO revision: the computations of
`agent_instructions_computations_2.md`, i.e. finer discretisations, smooth
worst-case measures and the figures of the talk. **This release changes published
numbers**: the certified optimality gap of the manuscript becomes sharper, see the
first item below. It is a minor and not a patch release for exactly that reason,
the output of the same command differs from that of 0.2.0.

### Changed
- `Instanz.t_bar(i)` is now the diameter of the numerical support of the envelope
  of species `i`, i.e. of `T^s`, instead of the diameter of the whole time
  window. The bound of Theorem "inner_convergence" is evaluated per species, as
  the theorem states it, and becomes sharper: `sum_s kappa^s_N` at
  `delta_N = 1e-4` drops from 7.67 to 5.09 and the certified gap from 14.4 % to
  12.0 %. A cell counts as part of the support if its envelope mass is at least
  `1e-12`. On the shipped data that threshold drops nothing, since
  `aggregate_matrix` sets densities below `1e-5` to exactly zero; the support, and
  with it `T_bar_s`, is identical on all refinements of the grid. For data without
  such a cutoff the threshold would disregard at most `8e-9` of probability mass,
  i.e. change the bound by less than `1e-6`.
- `model_2.solve_dro_model` reports the solution on the grid of the instance,
  which differs from the grid it was called with when the grid is refined.

### Added
- Grid refinement by a zero-order hold: `hilfsfunktionen.refine_matrix_index`,
  `refine_matrix`, `refine_vector`, the parameters `refinement_factor` and
  `refine_recompute`, and the command line options `--refinement_factor`,
  `--refine_recompute`. `delta_N = 1e-4 / r` for the default data set.
  With `r = 1` the refined path reproduces the data grid instance bit-identically.
  The moment bounds and the envelope are transferred from the data grid rather
  than recomputed, which keeps the ambiguity set unchanged; `--refine_recompute`
  recomputes them and exists only to quantify the difference.
- `check_refinement.py`: verifies that the envelope masses of the refined
  sub-cells sum to the mass of the data cell (relative deviation below `1e-15`)
  and that `r = 1` is the identity.
- `inner_lp.smoothest_worst_case`: two-stage selection of a worst-case measure of
  least total variation or of least curvature among all optimal ones, with the
  value of the adversarial problem fixed to its optimum. The vertex that the
  simplex method returns oscillates between zero and the envelope from cell to
  cell; that oscillation is an artefact of vertex selection.
- `slide_plot.py`: the chromatogram figure in the colours of the talk, one
  vector PDF per fractionation window, including the worst-case densities, plus
  a `--from_json` path that redraws a panel from stored data without Gurobi.
- `experiments.py r1`: the sweep over the refinement factors, writing one row at
  a time so that a partial result survives.
