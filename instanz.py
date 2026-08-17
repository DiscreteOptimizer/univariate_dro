"""Instance pre-processing shared by the MIP model and the analysis tools.

This module collects the derived quantities of an instance (particle masses,
normalised densities, moment bounds, variance bound, envelope) in one place.
The arithmetic is exactly the one that was inlined in ``model_2.solve_dro_model``
up to version 0.1.1; the model now calls :func:`baue_instanz` instead, so that
the auxiliary scripts (error bound, reference enumeration, worst-case plot)
work on bit-identical data.
"""

# general imports
from dataclasses import dataclass, field
from typing import List

import numpy as np

# auxiliary functions
import hilfsfunktionen as aux
# params
from params import Params

# scaling factors to increase numeric stability (as in model_2 up to v0.1.1)
FACTOR = 1e06
Q_FACTOR = 100.


@dataclass
class Instanz:
    """Derived data of one instance.

    Attributes
    ----------
    time_points
        The discretisation grid :math:`T_N` (including the right end point of
        the last cell), as read from the data and aggregated.
    zeit_diskret
        The grid width :math:`\\delta_N`.
    anzahl_prozess
        Number of grid cells, i.e. ``len(time_points) - 1``.
    q0
        Particle masses, scaled by :data:`Q_FACTOR`.
    matrix_nom, matrix_min, matrix_max
        Densities normalised to unit mass (probability densities).
    ret_time, ret_time_minus, ret_time_plus
        :math:`\\mu_s`, :math:`\\mu_{s,-}`, :math:`\\mu_{s,+}`.
    dict_var_nom, dict_var_min, dict_var_max
        Variances of the nominal, fast and slow density.
    schwank_var_global
        Largest relative variance deviation, the factor of the variance bound.
    sigma2_plus
        The variance bound :math:`\\sigma^{2,s}_+`.
    varianz_schranke
        Right-hand side of the second moment constraint as used in the model,
        i.e. :math:`\\sigma^{2,s}_+ - \\mu_+\\mu_- + \\delta_N^2/4`.
    schlauch_rtd
        The envelope :math:`\\bar\\rho^s` on the grid, as a density.
    a_s
        The coefficients :math:`(\\mathbb{1}_{S_w}(s) - R)` of the purity
        constraint, *without* the particle mass ``q0``.
    """

    time_points: List[float]
    zeit_diskret: float
    anzahl_prozess: int
    groessen: List[int]
    q0: np.ndarray
    matrix_nom: np.ndarray
    matrix_min: np.ndarray
    matrix_max: np.ndarray
    totm_desired: float
    ret_time: List[float]
    ret_time_minus: List[float]
    ret_time_plus: List[float]
    dict_var_nom: List[float]
    dict_var_min: List[float]
    dict_var_max: List[float]
    schwank_var_global: float
    sigma2_plus: List[float]
    varianz_schranke: List[float]
    schlauch_rtd: List[List[float]]
    a_s: List[float]
    factor: float = FACTOR
    q_factor: float = Q_FACTOR

    def a_coeff(self, i: int) -> float:
        """Coefficient :math:`a^s` of species ``i`` as it enters the model."""
        return self.a_s[i] * self.q0[i]

    def rho_max(self, i: int) -> float:
        """Maximum :math:`\\rho^{max,s}` of the envelope of species ``i``."""
        return max(self.schlauch_rtd[i])

    def t_bar(self) -> float:
        """:math:`\\bar T = \\max T - \\min T`."""
        return self.time_points[-1] - self.time_points[0]


def baue_instanz(time_points, matrix_nom_roh, matrix_min_roh, matrix_max_roh,
                 params: Params) -> Instanz:
    """Compute all derived quantities of an instance.

    The order of operations is identical to the one used in
    ``model_2.solve_dro_model`` up to version 0.1.1, so that results do not
    change.
    """
    # initialize particle mass vector
    q0 = np.zeros(len(matrix_nom_roh))
    # calculate mass of particle i
    for i in range(len(matrix_nom_roh)):
        q0[i] = aux.flaeche(time_points, matrix_nom_roh[i])

    # initialize scaled data
    matrix_nom = np.zeros((len(matrix_nom_roh), len(matrix_nom_roh[0])))
    matrix_min = np.zeros((len(matrix_nom_roh), len(matrix_nom_roh[0])))
    matrix_max = np.zeros((len(matrix_nom_roh), len(matrix_nom_roh[0])))

    # scale data...
    for i, (el1, el2, el3, el4) in enumerate(zip(q0, matrix_nom_roh, matrix_min_roh, matrix_max_roh)):
        matrix_nom[i] = el2 / el1
        matrix_min[i] = el3 / el1
        matrix_max[i] = el4 / el1

    # particle masses
    for i in range(len(q0)):
        q0[i] *= Q_FACTOR

    # particle indices
    groessen = list(i for i in range(len(matrix_nom)))
    # number of time steps
    anzahl_prozess = len(time_points) - 1

    # total time interval length
    aux_sum = 0
    for i in range(len(time_points) - 1):
        aux_sum += time_points[i + 1] - time_points[i]

    # calculate \delta_N as the mean of all time steps
    zeit_diskret = aux_sum / (len(time_points) - 1)

    # total mass of desired peak is 1...
    totm_desired = aux.flaeche(time_points, matrix_nom[params.wunschgroesse])

    # mu, mu^-, mu^+
    ret_time = aux.baue_mu_list(time_points, matrix_nom)
    ret_time_minus = aux.baue_mu_list(time_points, matrix_min)
    ret_time_plus = aux.baue_mu_list(time_points, matrix_max)

    # some intermediate results for the variance bound
    dict_var_nom = aux.baue_var_list(time_points, matrix_nom, ret_time)
    dict_var_min = aux.baue_var_list(time_points, matrix_min, ret_time_minus)
    dict_var_max = aux.baue_var_list(time_points, matrix_max, ret_time_plus)

    # some intermediate results for the variance bound
    schwankung_min = [abs(mini / nomi) for mini, nomi in zip(dict_var_min, dict_var_nom)]
    schwankung_max = [abs(maxi / nomi) for maxi, nomi in zip(dict_var_max, dict_var_nom)]
    schwank_var_global = max([max(schwankung_min), max(schwankung_max)])

    # variance bound \sigma^{2,s}_+
    sigma2_plus = [schwank_var_global * nomi for nomi in dict_var_nom]

    # right-hand side of the second moment constraint of the model
    varianz_schranke = [schwank_var_global * nomi - muplus * muminus + zeit_diskret ** 2 / 4.
                        for nomi, muplus, muminus in zip(dict_var_nom, ret_time_plus, ret_time_minus)]

    # calculate envelopes
    schlauch_rtd = aux.schlauch_chromatogramm(matrix_nom, matrix_min, matrix_max)

    # optionally replace the rectangle rule for the envelope masses by the exact
    # integral over the cell (verification task V2). The model multiplies
    # schlauch_rtd by zeit_diskret, so the exact mass is divided by it here.
    if getattr(params, 'exact_envelope_mass', False):
        import envelope as env_mod
        for i in groessen:
            exact = env_mod.exact_cell_masses(time_points,
                                             [matrix_nom[i], matrix_min[i], matrix_max[i]])
            for j in range(anzahl_prozess):
                schlauch_rtd[i][j] = float(exact[j]) / zeit_diskret

    # purity coefficients (without particle mass)
    a_s = []
    for i in groessen:
        if i == params.wunschgroesse:
            a_s.append(1 - params.reinheit)
        else:
            a_s.append(-params.reinheit)

    return Instanz(
        time_points=time_points,
        zeit_diskret=zeit_diskret,
        anzahl_prozess=anzahl_prozess,
        groessen=groessen,
        q0=q0,
        matrix_nom=matrix_nom,
        matrix_min=matrix_min,
        matrix_max=matrix_max,
        totm_desired=totm_desired,
        ret_time=ret_time,
        ret_time_minus=ret_time_minus,
        ret_time_plus=ret_time_plus,
        dict_var_nom=dict_var_nom,
        dict_var_min=dict_var_min,
        dict_var_max=dict_var_max,
        schwank_var_global=schwank_var_global,
        sigma2_plus=sigma2_plus,
        varianz_schranke=varianz_schranke,
        schlauch_rtd=schlauch_rtd,
        a_s=a_s,
    )
