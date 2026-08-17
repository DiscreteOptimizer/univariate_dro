"""Regeneration of the residence time distributions for other values of eps_ACN.

The shipped data sets only cover three values of the acetonitrile uncertainty,

    daten/Small_E  : eps_ACN = 0.0040   (grid 1e-3)
    daten/Medium_E : eps_ACN = 0.0042   (grid 1e-3)
    daten/Large_E  : eps_ACN = 0.0044   (grid 1e-3)
    daten/Long_E   : eps_ACN = 0.0042   (grid 1e-4, same distributions as Medium_E)

so that the sweep over eps_ACN of the manuscript cannot be produced from them.
This module regenerates the data for an arbitrary eps_ACN. It rests on two
observations about the shipped files, both of which are verified by
``--validate``:

1. Every shipped density is exactly a normal density with mean ``mu`` and
   standard deviation ``mu / sqrt(NTP)``, scaled by the particle mass ``q0_s``,
   evaluated on the grid and truncated to ``T = [2.75, 3.7]`` (the truncation
   removes less than 1e-16 of the mass). This is the model of the manuscript,
   Table "parameter" and Example "envelope_normal".

2. The only ingredient that is *not* contained in the code is the map
   ``r_ACN -> mu_s(r_ACN)`` of Supper et al. It is recovered from the shipped
   data: the four data sets provide the seven points
   ``r_ACN in {0.25} u {0.25 +- eps : eps in {0.004, 0.0042, 0.0044}}``
   per species, and ``log mu_s`` is interpolated through them by a cubic
   polynomial. The interpolation residual on the seven points is below 3e-9
   minutes; the difference between the cubic and the quadratic interpolant at
   eps_ACN = 0.0063, i.e. outside the calibration range, is below 1e-5 minutes.
   Numbers obtained for eps_ACN outside [0.2456, 0.2544] therefore rest on an
   extrapolation of a surrogate of the model of Supper et al., not on that model
   itself, and must be labelled as such.

Usage::

    python3.12 generate_data.py --validate
    python3.12 generate_data.py --eps 0.0021 0.0063        # writes daten/Gen_0.0021_E, ...
    python3.12 run_funktionen.py --sample Gen_0.0021_E

The written directories follow the naming convention of the shipped data, so
that ``--sample Gen_<eps>_E`` picks them up.
"""

import argparse
import glob
import os

import numpy as np

# number of theoretical plates, Table "parameter" of the manuscript
NTP = 120000.
# nominal acetonitrile ratio
R_NOM = 0.25
# eps_ACN of the shipped data sets, established by comparing the resulting
# moment bounds with Table "inputvalues" of the manuscript and with the
# commented-out table in the source of comp-results.tex
SHIPPED_EPS = {'Small_E': 0.0040, 'Medium_E': 0.0042, 'Large_E': 0.0044, 'Long_E': 0.0042}
# the data set the grid and the particle masses are taken from
REFERENCE_SET = 'Long_E'
SPECIES = (30, 31, 32, 33)


def _read(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]


def _moments(t, y):
    mass = float(np.sum(y[:-1] * np.diff(t)))
    mu = float(np.sum(t * y) / np.sum(y))
    return mass, mu


def read_calibration(daten='daten'):
    """Recover ``mu_s(r_ACN)`` and the particle masses from the shipped data."""
    points = {s: {} for s in SPECIES}
    for name, eps in SHIPPED_EPS.items():
        if name == REFERENCE_SET:
            continue                      # duplicate of Medium_E, no new information
        for s in SPECIES:
            for kind, r in (('', R_NOM), ('_min', R_NOM + eps), ('_max', R_NOM - eps)):
                t, y = _read(os.path.join(daten, name + kind, f'{name}{kind}_{s}.txt'))
                points[s][round(r, 6)] = _moments(t, y)[1]

    fits = {}
    for s in SPECIES:
        r = np.array(sorted(points[s]))
        mu = np.array([points[s][x] for x in r])
        cubic = np.polyfit(r, np.log(mu), 3)
        quadratic = np.polyfit(r, np.log(mu), 2)
        resid = float(np.abs(np.exp(np.polyval(cubic, r)) - mu).max())
        fits[s] = {'cubic': cubic, 'quadratic': quadratic, 'residual': resid,
                   'r': r, 'mu': mu}

    # grid and particle masses from the reference data set
    grid = None
    masses = {}
    amp_resid = {}
    for s in SPECIES:
        t, y = _read(os.path.join(daten, REFERENCE_SET, f'{REFERENCE_SET}_{s}.txt'))
        grid = t if grid is None else grid
        mass, mu = _moments(t, y)
        sd = mu / np.sqrt(NTP)
        pdf = np.exp(-0.5 * ((t - mu) / sd) ** 2) / (sd * np.sqrt(2. * np.pi))
        sel = pdf > 1e-6 * pdf.max()
        amp = y[sel] / pdf[sel]
        masses[s] = float(np.median(amp))
        amp_resid[s] = float(np.abs(amp / masses[s] - 1.).max())
    return {'fits': fits, 'grid': grid, 'masses': masses, 'amplitude_residual': amp_resid}


def mu_of_r(cal, s, r, kind='cubic'):
    """Interpolated mean retention time of species ``s`` at ratio ``r``."""
    return float(np.exp(np.polyval(cal['fits'][s][kind], r)))


def density(t, mu, mass):
    """Scaled normal density with standard deviation ``mu / sqrt(NTP)``."""
    sd = mu / np.sqrt(NTP)
    return mass * np.exp(-0.5 * ((t - mu) / sd) ** 2) / (sd * np.sqrt(2. * np.pi))


def write_instance(cal, eps, name=None, daten='daten', grid=None, kind='cubic'):
    """Write the three directories of one instance and return their names."""
    name = name if name is not None else f'Gen_{eps:g}_E'
    grid = cal['grid'] if grid is None else grid
    written = []
    info = {'eps_ACN': eps, 'name': name, 'interpolation': kind, 'species': {}}
    for suffix, r in (('', R_NOM), ('_min', R_NOM + eps), ('_max', R_NOM - eps)):
        d = os.path.join(daten, name + suffix)
        os.makedirs(d, exist_ok=True)
        for s in SPECIES:
            mu = mu_of_r(cal, s, r, kind)
            y = density(grid, mu, cal['masses'][s])
            path = os.path.join(d, f'{name}{suffix}_{s}.txt')
            with open(path, 'w') as f:
                for ti, yi in zip(grid, y):
                    f.write(f"{float(ti)!r} {float(yi)!r}\n")
            info['species'].setdefault(s, {})[suffix or 'nom'] = mu
        written.append(d)
    return written, info


def validate(cal, daten='daten', tol=1e-6):
    """Regenerate every shipped data set and compare with the shipped files.

    The deviation is dominated by the residual of the interpolation of
    ``mu_s(r_ACN)``: an error of 1e-8 minutes in the mean moves the density on the
    steep flank of the peak by about 1e-7 of its maximum.
    """
    print(f"calibration: interpolation residuals "
          f"{ {s: '%.2e' % cal['fits'][s]['residual'] for s in SPECIES} }")
    print(f"             amplitude residuals   "
          f"{ {s: '%.2e' % cal['amplitude_residual'][s] for s in SPECIES} }")
    worst = 0.
    worst_mu = 0.
    for name, eps in SHIPPED_EPS.items():
        for suffix, r in (('', R_NOM), ('_min', R_NOM + eps), ('_max', R_NOM - eps)):
            for s in SPECIES:
                path = os.path.join(daten, name + suffix, f'{name}{suffix}_{s}.txt')
                t, y = _read(path)
                mu = mu_of_r(cal, s, r)
                y_gen = density(t, mu, cal['masses'][s])
                scale = max(y.max(), 1e-30)
                err = float(np.abs(y - y_gen).max() / scale)
                d_mu = abs(mu - _moments(t, y)[1])
                worst = max(worst, err)
                worst_mu = max(worst_mu, d_mu)
                print(f"  {name+suffix:>16} s={s}: max rel. deviation {err:.2e}, "
                      f"deviation of the mean {d_mu:.2e} min")
    print(f"worst relative deviation over all shipped files: {worst:.3e} "
          f"({'ok' if worst < tol else 'CHECK'})")
    print(f"worst deviation of the mean retention time: {worst_mu:.3e} minutes")
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--daten', default='daten')
    parser.add_argument('--eps', type=float, nargs='*', default=[],
                        help="Values of eps_ACN to generate.")
    parser.add_argument('--name', default=None,
                        help="Directory name for a single --eps value.")
    parser.add_argument('--validate', action='store_true',
                        help="Regenerate the shipped data sets and report the deviation.")
    parser.add_argument('--table', action='store_true',
                        help="Print the resulting moment bounds without writing files.")
    args = parser.parse_args()

    cal = read_calibration(args.daten)
    if args.validate:
        validate(cal, args.daten)
    if args.table:
        for eps in args.eps:
            print(f"eps_ACN = {eps:g}")
            for s in SPECIES:
                lo = mu_of_r(cal, s, R_NOM + eps)
                hi = mu_of_r(cal, s, R_NOM - eps)
                lo2 = mu_of_r(cal, s, R_NOM + eps, 'quadratic')
                hi2 = mu_of_r(cal, s, R_NOM - eps, 'quadratic')
                print(f"  s={s}: mu_-={lo:.6f} mu_+={hi:.6f} half width={(hi-lo)/2:.6f} "
                      f"(cubic-quadratic: {abs(lo-lo2):.1e}/{abs(hi-hi2):.1e})")
    if args.eps and not args.table:
        for eps in args.eps:
            written, info = write_instance(cal, eps, args.name if len(args.eps) == 1 else None,
                                          args.daten)
            print(f"eps_ACN = {eps:g} -> {', '.join(written)}")
            for s, mus in info['species'].items():
                print(f"    s={s}: mu={mus['nom']:.6f} mu_-={mus['_min']:.6f} "
                      f"mu_+={mus['_max']:.6f}")


if __name__ == "__main__":
    main()
