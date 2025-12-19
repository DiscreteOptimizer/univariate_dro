# general imports
from dataclasses import dataclass

# ----------------------------
# Parameter-Dataclass
# ----------------------------

@dataclass
class Params:
    aggregation_factor: int
    reinheit: float
    wunschgroesse: int
    sample: str
    fix: bool
    fix_lower: int
    fix_upper: int
    nominal: bool

