"""Exact envelope masses versus the rectangle rule (verification task V2).

The discretised envelope constraint of the manuscript is

    P([tau, tau + delta_N]) <= rho_bar_N^+(tau) := int_tau^{tau+delta_N} rho_bar(t) dt

for ``eta_N = 0``, whereas the implementation uses the rectangle rule
``rho_bar_N^+(tau) = delta_N * rho_bar(tau)`` with the value at the *left* end
point of the cell. On a cell where the envelope increases, the rectangle rule
*underestimates* the true envelope mass, so the implemented constraint is
stronger than the one of the manuscript and excludes measures that belong to the
ambiguity set (Eq. "ambiguity_set"). On a cell where the envelope decreases it is
the other way round. The two ambiguity sets are therefore not comparable, and in
particular the implemented one is not a relaxation of the ambiguity set of the
manuscript.

This module computes the exact cell masses, so that the size of the effect can be
quantified and the model can be re-solved with the correct masses. It exploits
that the shipped densities are exactly scaled normal densities with standard
deviation ``mu / sqrt(NTP)``, see ``generate_data.py``; ``verify_fit`` reports the
deviation between the fitted analytic densities and the data.

The envelope itself is the one of ``hilfsfunktionen.schlauch_chromatogramm``:
outside the interval spanned by the three peak positions it is the pointwise
maximum of the three densities, inside it is constant and equal to the largest of
the three peak values.

Only the native data grid is supported (``aggregation_factor = 1``); on an
aggregated grid the entries are cell averages and no longer point values of a
normal density.
"""

import numpy as np

# number of theoretical plates, Table "parameter" of the manuscript
NTP = 120000.

# 16-point Gauss-Legendre rule on [-1, 1]
_GL_X, _GL_W = np.polynomial.legendre.leggauss(16)


def _fit(t, y):
    """Fit a scaled normal density with sd = mu / sqrt(NTP) to the samples."""
    mu = float(np.sum(t * y) / np.sum(y))
    sd = mu / np.sqrt(NTP)
    pdf = np.exp(-0.5 * ((t - mu) / sd) ** 2) / (sd * np.sqrt(2. * np.pi))
    sel = pdf > 1e-8 * pdf.max()
    scale = float(np.median(y[sel] / pdf[sel]))
    resid = float(np.abs(y[sel] - scale * pdf[sel]).max() / y.max())
    return {'mu': mu, 'sd': sd, 'scale': scale, 'residual': resid}


def fit_species(time_points, rows):
    """Analytic description of the three densities (nominal, fast, slow)."""
    t = np.asarray(time_points, dtype=float)
    return [_fit(t, np.asarray(r, dtype=float)) for r in rows]


def _density(fits, t):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.zeros((len(fits), t.size))
    for k, f in enumerate(fits):
        out[k] = f['scale'] * np.exp(-0.5 * ((t - f['mu']) / f['sd']) ** 2) \
            / (f['sd'] * np.sqrt(2. * np.pi))
    return out


def envelope_description(time_points, rows):
    """Fits, flat top value and the interval on which the envelope is constant.

    ``rows`` are the nominal, fast and slow density of one species, sampled on
    ``time_points``.
    """
    t = np.asarray(time_points, dtype=float)
    rows = [np.asarray(r, dtype=float) for r in rows]
    fits = fit_species(t, rows)
    args = [int(np.argmax(r)) for r in rows]
    lo, hi = min(args), max(args)
    flat = max(float(r.max()) for r in rows)
    return {'fits': fits, 'flat': flat, 'flat_from': float(t[lo]), 'flat_to': float(t[hi]),
            'index_from': lo, 'index_to': hi}


def analytic_envelope(desc, t):
    """The envelope of ``desc`` evaluated at ``t``."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    val = _density(desc['fits'], t).max(axis=0)
    inside = (t >= desc['flat_from']) & (t <= desc['flat_to'])
    val[inside] = desc['flat']
    return val


def verify_fit(time_points, rows, grid_envelope):
    """Largest deviation between the analytic envelope and the grid envelope."""
    t = np.asarray(time_points, dtype=float)
    grid = np.asarray(grid_envelope, dtype=float)
    ana = analytic_envelope(envelope_description(t, rows), t)
    scale = max(grid.max(), 1e-30)
    return float(np.abs(ana - grid).max() / scale)


def _integrate(desc, a, b):
    """Integral of the envelope over ``[a, b]``, split at the two kinks."""
    nodes = sorted({a, b} | {x for x in (desc['flat_from'], desc['flat_to']) if a < x < b})
    total = 0.
    for lo, hi in zip(nodes[:-1], nodes[1:]):
        if hi <= lo:
            continue
        mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
        total += half * float(np.dot(_GL_W, analytic_envelope(desc, mid + half * _GL_X)))
    return total


def exact_cell_masses(time_points, rows):
    """``int_{t_j}^{t_j+delta_N} rho_bar`` for every grid cell of one species."""
    t = np.asarray(time_points, dtype=float)
    desc = envelope_description(t, rows)
    return np.array([_integrate(desc, t[j], t[j + 1]) for j in range(len(t) - 1)])


def species_rows(inst, i):
    """The three densities of species ``i`` of an instance."""
    return [inst.matrix_nom[i], inst.matrix_min[i], inst.matrix_max[i]]


def exact_cell_masses_inst(inst, i):
    return exact_cell_masses(inst.time_points, species_rows(inst, i))


def comparison(inst):
    """Rectangle rule versus exact envelope masses, per species."""
    d = inst.zeit_diskret
    rows = []
    for i in inst.groessen:
        exact = exact_cell_masses_inst(inst, i)
        rect = d * np.asarray(inst.schlauch_rtd[i][:inst.anzahl_prozess])
        diff = exact - rect
        rows.append({
            'species_index': i,
            'fit_residual': verify_fit(inst.time_points, species_rows(inst, i),
                                       inst.schlauch_rtd[i]),
            'sum_exact': float(exact.sum()),
            'sum_rectangle': float(rect.sum()),
            'max_positive_diff': float(diff.max()),
            'max_negative_diff': float(diff.min()),
            'sum_positive_diff': float(diff[diff > 0].sum()),
            'sum_negative_diff': float(diff[diff < 0].sum()),
            'max_relative_diff': float(np.abs(diff).max() / rect.max()),
            'n_cells_rectangle_too_small': int((diff > 0).sum()),
            'n_cells_rectangle_too_large': int((diff < 0).sum()),
        })
    return rows
