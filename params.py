# general imports
from dataclasses import dataclass
from typing import Optional

# ----------------------------
# Parameter-Dataclass
# ----------------------------

@dataclass
class Params:
    """Parameters of one run.

    The first block is the interface of version 0.1.1 and unchanged. The second
    block was added in version 0.2.0 for the additional experiments; all of its
    entries have defaults that reproduce the behaviour of version 0.1.1.
    """

    aggregation_factor: int
    reinheit: float
    wunschgroesse: int
    sample: str
    fix: bool
    fix_lower: int
    fix_upper: int
    nominal: bool

    # --- added in 0.2.0 -----------------------------------------------------
    # right-hand side of the purity constraint (32b). 0.0 is the safe
    # approximation; a negative value -sum_s Delta_N^s yields the upper bound
    # val_N^up of Corollary "enclosure".
    purity_rhs: float = 0.0
    # drop the relaxed second moment constraint by fixing the dual variable
    # y_5 = dualvariablen[i][4] to zero
    no_second_moment: bool = False
    # write the Gurobi log to this file in addition to the console
    log_file: Optional[str] = None
    # show/save the matplotlib figure of run_funktionen.run
    plot: bool = True
    # write the solver statistics and the solution as JSON to this file
    stats_file: Optional[str] = None
    # keep the model object and the variable handles in the result dictionary
    keep_model: bool = False
    # use the exact envelope mass of every grid cell instead of the rectangle
    # rule delta_N * rho_bar(tau) used up to version 0.1.1; see envelope.py.
    # Only available on the native data grid (aggregation_factor = 1).
    exact_envelope_mass: bool = False
    # optimise only the length of the fractionation interval instead of the
    # hierarchy of three objectives. The optimal length is the same, but a single
    # objective gives access to MIPGap and ObjBound, which is needed to report a
    # valid bound when a run is stopped by the time limit.
    single_objective: bool = False
    # Gurobi time limit in seconds; None means no limit
    time_limit: Optional[float] = None
