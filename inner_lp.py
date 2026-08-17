"""Linear programs for the adversarial problem of one species.

Two linear programs are provided for the inner problem

    val^s(x^-, x^+) = min_{P in U^s} a^s P([x^-, x^+])                     (25)

on a fixed grid ``T_N``:

``LowerLP``
    Problem "safe_approx_primal", i.e. the re-dualisation of the discretised
    dual problem. Its optimal value is a *lower* bound on ``val^s``; this is the
    relaxation on which the safe approximation of the manuscript is built. It
    uses the grid-rounded indicator ``1^c`` of equations "taudef"/"taudef2",
    i.e. one cell is removed at each end for ``a^s > 0`` and one cell is added at
    each end for ``a^s < 0``.

``UpperLP``
    A *restriction* of (25) to measures whose density is constant on every grid
    cell and bounded by the envelope on the whole cell. Every feasible point of
    this LP is a feasible point of (25), hence its optimal value is an *upper*
    bound on ``val^s``. First and second moment of such a measure are evaluated
    exactly, and the objective uses the exact interval ``[x^-, x^+]``.

Together the two LPs bracket the true value ``val^s``. The bracket is much
sharper than the a priori bound ``Delta_N^s`` of ``error_bound.py``, and it is
what makes the reference solution of ``reference.py`` a certified enclosure
rather than a point estimate.

Both classes keep one Gurobi model alive and only update the objective
coefficients that change when the fractionation interval moves, so that a scan
over many intervals is solved by a few dual simplex iterations per interval.

The interval is always given by two grid indices ``(lo, hi)`` with
``x^- = time_points[lo]``, ``x^+ = time_points[hi]``.
"""

import gurobipy as gp
import numpy as np

from instanz import Instanz


class _BaseLP:
    def __init__(self, inst: Instanz, i: int, env: gp.Env = None):
        self.inst = inst
        self.i = i
        self.a = inst.a_coeff(i)
        self.delta = inst.zeit_diskret
        self.t = np.asarray(inst.time_points, dtype=float)
        self.n = len(self.t)
        self.K = inst.anzahl_prozess          # number of grid cells
        self.mu_minus = inst.ret_time_minus[i]
        self.mu_plus = inst.ret_time_plus[i]
        self.sigma2_plus = inst.sigma2_plus[i]
        self.env = np.asarray(inst.schlauch_rtd[i], dtype=float)
        self.m = gp.Model(f"inner_{type(self).__name__}_{i}", env=env) if env is not None \
            else gp.Model(f"inner_{type(self).__name__}_{i}")
        self.m.setParam('OutputFlag', 0)
        self.m.setParam('FeasibilityTol', 1e-9)
        self.m.setParam('OptimalityTol', 1e-9)
        self._current = None                  # currently set index set
        self._build()

    def _build(self):
        raise NotImplementedError

    def _coeff_vector(self, lo, hi):
        raise NotImplementedError

    def solve(self, lo, hi):
        """Optimal value of the LP for the interval ``[t[lo], t[hi]]``."""
        self._set_objective(lo, hi)
        self.m.optimize()
        if self.m.status != gp.GRB.Status.OPTIMAL:
            raise RuntimeError(f"inner LP not solved to optimality: status {self.m.status}")
        return self.m.ObjVal

    @property
    def iter_count(self):
        return self.m.IterCount


class LowerLP(_BaseLP):
    """Lower bound on ``val^s`` (Problem "safe_approx_primal")."""

    def _build(self):
        m, K, t, d = self.m, self.K, self.t, self.delta
        # w^-_j sits at t[j] for j = 0..K-1, w^+_j sits at t[j+1] for j = 0..K-2
        self.wm = m.addVars(K, lb=0., ub=gp.GRB.INFINITY, name='wm')
        self.wp = m.addVars(K - 1, lb=0., ub=gp.GRB.INFINITY, name='wp')

        zweimu = self.mu_minus + self.mu_plus
        pos_m = t[:K]
        pos_p = t[1:K]

        mass = gp.quicksum(self.wm[j] for j in range(K)) + gp.quicksum(self.wp[j] for j in range(K - 1))
        m.addConstr(mass == 1., name='mass')

        mean = gp.quicksum(pos_m[j] * self.wm[j] for j in range(K)) \
            + gp.quicksum(pos_p[j] * self.wp[j] for j in range(K - 1))
        m.addConstr(mean >= self.mu_minus, name='mean_lo')
        m.addConstr(mean <= self.mu_plus, name='mean_hi')

        # second moment, cf. constraint "variance_ref"
        var = gp.quicksum((pos_m[j] ** 2 - zweimu * pos_m[j]) * self.wm[j] for j in range(K)) \
            + gp.quicksum((pos_p[j] ** 2 - zweimu * pos_p[j]) * self.wp[j] for j in range(K - 1))
        rhs = self.sigma2_plus - self.mu_plus * self.mu_minus + d ** 2 / 4.
        m.addConstr(var <= rhs, name='variance')

        # envelope mass of cell j, evaluated by the rectangle rule (eta_N = 0)
        for j in range(K):
            terms = self.wm[j] + (self.wp[j] if j < K - 1 else 0.)
            m.addConstr(terms <= d * self.env[j], name=f'env_{j}')
        m.update()

    # Number of cells by which the interval is shrunk at its left and right end
    # (a > 0) resp. widened (a < 0). The default (1, 1) is the grid-rounded
    # indicator 1^c of equations "taudef"/"taudef2", where the widening for
    # a < 0 is by one more cell on the right than the manuscript prescribes, so
    # that the resulting value is a lower bound on the value of the manuscript's
    # linear program as well.
    #
    # The MIP of model_2.py uses the further weakened indicator of Systems
    # (29)/(30). There, b_tilde equals one on the grid indices [lo, hi-1], hence
    # fract_p (a > 0) is one on [lo+2, hi-3] and fract_n (a < 0) on [lo-1, hi].
    # Those are the values used with mip_indicator, and they reproduce the inner
    # values of the MIP.
    shift_pos = (1, 1)
    shift_neg = (1, 1)

    def support_indices(self, lo, hi):
        """Indices of the cells whose weights carry the objective coefficient."""
        if self.a > 0:
            # 1^c is zero outside [x^-, x^+] and one on [tau^-_N, tau^+_N]
            first, last = lo + self.shift_pos[0], hi - self.shift_pos[1]
        else:
            # 1^c is one on [x^-, x^+] and zero outside [tau^-_N, tau^+_N]
            first, last = lo - self.shift_neg[0], hi + self.shift_neg[1]
        first = max(first, 0)
        last = min(last, self.K - 1)
        if first > last:
            return range(0, 0)
        return range(first, last + 1)

    def _set_objective(self, lo, hi):
        new = set(self.support_indices(lo, hi))
        old = self._current if self._current is not None else set()
        for j in new - old:
            self.wm[j].Obj = self.a
            if j < self.K - 1:
                self.wp[j].Obj = self.a
        for j in old - new:
            self.wm[j].Obj = 0.
            if j < self.K - 1:
                self.wp[j].Obj = 0.
        self._current = new


class UpperLP(_BaseLP):
    """Upper bound on ``val^s``: restriction to cell-wise constant densities."""

    def _build(self):
        m, K, t, d = self.m, self.K, self.t, self.delta
        # density level of cell j is p[j]/delta and must not exceed the envelope
        # anywhere in the cell. The envelope is piecewise monotone or has an
        # interior maximum on a cell, so its minimum over the cell is attained at
        # an end point.
        cap = d * np.minimum(self.env[:K], self.env[1:K + 1])
        self.p = m.addVars(K, lb=0., ub=[float(c) for c in cap], name='p')

        mid = t[:K] + d / 2.
        zweimu = self.mu_minus + self.mu_plus

        m.addConstr(gp.quicksum(self.p[j] for j in range(K)) == 1., name='mass')
        mean = gp.quicksum(mid[j] * self.p[j] for j in range(K))
        m.addConstr(mean >= self.mu_minus, name='mean_lo')
        m.addConstr(mean <= self.mu_plus, name='mean_hi')
        # exact second moment of a uniform density on a cell: mid^2 + delta^2/12
        var = gp.quicksum((mid[j] ** 2 + d ** 2 / 12. - zweimu * mid[j]) * self.p[j] for j in range(K))
        m.addConstr(var <= self.sigma2_plus - self.mu_plus * self.mu_minus, name='variance')
        m.update()

    # If ``cover`` is set, the objective covers every interval [y^-, y^+] with
    # |y^- - t[lo]| <= delta_N and |y^+ - t[hi]| <= delta_N: the desired species
    # (a > 0) then uses the cells of the largest such interval and the
    # contaminants (a < 0) the cells of the smallest one. The resulting value is
    # an upper bound on val^s(y^-, y^+) for every such pair, which makes a scan
    # over grid points a valid test for all real intervals in between.
    cover = False

    def support_indices(self, lo, hi):
        """Cells that are contained in the exact interval ``[x^-, x^+]``."""
        if self.cover:
            if self.a > 0:
                lo, hi = lo - 1, hi + 1
            else:
                lo, hi = lo + 1, hi - 1
        return range(max(lo, 0), min(max(hi, 0), self.K))

    def _set_objective(self, lo, hi):
        new = set(self.support_indices(lo, hi))
        old = self._current if self._current is not None else set()
        for j in new - old:
            self.p[j].Obj = self.a
        for j in old - new:
            self.p[j].Obj = 0.
        self._current = new

    def measure(self):
        """The optimal density levels of the last solve (worst-case measure)."""
        return np.array([self.p[j].X for j in range(self.K)]) / self.delta


def build_all(inst: Instanz, kind='lower', env=None, mip_indicator=False):
    """One LP per species.

    With ``mip_indicator`` the lower bound LPs use the weakened indicator of the
    mixed-integer model instead of ``1^c``, i.e. they reproduce the inner values
    of the MIP exactly.
    """
    cls = LowerLP if kind == 'lower' else UpperLP
    lps = [cls(inst, i, env=env) for i in inst.groessen]
    if mip_indicator and kind == 'lower':
        for lp in lps:
            lp.shift_pos, lp.shift_neg = (2, 3), (1, 0)
    return lps


def purity_from_values(inst: Instanz, values):
    """Worst-case purity that corresponds to the per-species values ``val^s``.

    ``val^s / a^s`` is the worst-case probability ``P^s([x^-,x^+])``: the minimum
    for the desired species (``a^s > 0``) and the maximum for the contaminants
    (``a^s < 0``). The worst-case purity is monotone in these quantities, cf. the
    derivation of Problem "DRO_chromatography" in the manuscript.
    """
    probs = [values[i] / inst.a_coeff(i) for i in inst.groessen]
    masses = [inst.q0[i] * probs[i] for i in inst.groessen]
    total = sum(masses)
    if total <= 0.:
        return None
    desired = masses[[i for i in inst.groessen
                      if inst.a_s[i] > 0][0]]
    return desired / total
