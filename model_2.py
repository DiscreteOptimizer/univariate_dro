# general imports
import json
import time

# gurobi
import gurobipy as gp

# auxiliary functions
import hilfsfunktionen as aux
# instance pre-processing
from instanz import baue_instanz, Instanz
# params
from params import Params


def build_dro_model(inst: Instanz, params: Params):
    """Build the MIP (Prob: MIP_onedim) for the instance ``inst``.

    Returns the Gurobi model and a dictionary of variable handles and objective
    expressions. The formulation is the one of version 0.1.1; the only additions
    are the two switches ``params.purity_rhs`` (right-hand side of Constraint
    (32b), used for the enclosure of Corollary "enclosure") and
    ``params.no_second_moment`` (fixes the dual variable of the relaxed second
    moment constraint to zero).
    """
    time_points = inst.time_points
    groessen = inst.groessen
    anzahl_prozess = inst.anzahl_prozess
    zeit_diskret = inst.zeit_diskret
    factor = inst.factor
    q0 = inst.q0
    matrix_nom = inst.matrix_nom
    ret_time_minus = inst.ret_time_minus
    ret_time_plus = inst.ret_time_plus
    varianz_schranke = inst.varianz_schranke
    schlauch_rtd = inst.schlauch_rtd
    a_s = inst.a_s

    m = gp.Model("chromatogram_dro")

    ##########################
    # VARIABLES ##############
    ##########################

    # \tilde{b} vars and modified \tilde{b}^- (for system (29)), \tilde{b}^+ (for system (31))
    fractionierung = []
    fract_p = []
    fract_n = []
    for i in range(anzahl_prozess + 1):
        fractionierung.append(m.addVar(vtype="B", name='fractionierung_%i' % i))
        fract_p.append(m.addVar(vtype="B", name='fract_p_%i' % i))
        fract_n.append(m.addVar(vtype="B", name='fract_n_%i' % i))

    # y variables
    dualvariablen = []
    for i in groessen:
        dualvariablen_zeile = []
        for j in range(5):
            dualvariablen_zeile.append(m.addVar(name='dualvar_%i_%i' % (i, j)))  # s. (20n)
        dualvariablen.append(dualvariablen_zeile)

    # z variables
    schlauchvariable = []
    for i in groessen:
        hilfe = []
        for j in range(anzahl_prozess):
            hilfe.append(m.addVar(name='envvar_%i_%i' % (i, j)))    # s. (20n)
        schlauchvariable.append(hilfe)

    # \Delta^-, \Delta^+ variables
    jump_variable = []
    for t in range(anzahl_prozess + 1):
        jump_variable.append([m.addVar(vtype="B", name='jump_pos_%i' % (t)), m.addVar(vtype="B", name='jump_neg_%i' % (t))])

    ##########################
    # CONSTRAINTS ############
    ##########################

    # Modified systems (29) and (31)
    # first step
    m.addConstr(fractionierung[0] == jump_variable[0][0] - jump_variable[0][1], name="fract_-1")
    # remaining steps
    for t in range(anzahl_prozess):
        m.addConstr(fractionierung[t + 1] == fractionierung[t] + jump_variable[t + 1][0] - jump_variable[t + 1][1], name=f"fract_{t}")

    # fract_p logic - force it to 0
    m.addConstr(fract_p[0] == 0, name="fract_p_0")
    m.addConstr(fract_p[1] == 0, name="fract_p_1")
    for t in range(2, anzahl_prozess-1):
        m.addConstr(fract_p[t] <= fractionierung[t - 2], name=f"fract_p_{t}a")
        m.addConstr(fract_p[t] <= fractionierung[t + 2], name=f"fract_p_{t}b")
    m.addConstr(fract_p[anzahl_prozess-1] == 0, name=f"fract_p_{anzahl_prozess-1}")
    m.addConstr(fract_p[anzahl_prozess] == 0, name=f"fract_p_{anzahl_prozess}")

    # fract_n logic - force it to 1
    for t in range(anzahl_prozess):
        m.addConstr(fract_n[t] >= fractionierung[t+1], name=f"fract_n_{t}a")
    for t in range(1, anzahl_prozess+1):
        m.addConstr(fract_n[t] >= fractionierung[t - 1], name=f"fract_n_{t}b")

    # (27a), (27b)
    m.addConstr(gp.quicksum(jump_variable[t][0] for t in range(anzahl_prozess + 1)) == 1, name="jump_0")
    m.addConstr(gp.quicksum(jump_variable[t][1] for t in range(anzahl_prozess + 1)) == 1, name="jump_1")

    # fix fractionation times if params.fix
    if params.fix:
        m.addConstr(jump_variable[params.fix_lower][0] == 1., name="jump_0_fix")
        m.addConstr(jump_variable[params.fix_upper][1] == 1., name="jump_0_fix")

    # drop the relaxed second moment constraint (12)/(Eq: Sec2_second_moment_true)
    # by fixing its dual variable to zero
    if params.no_second_moment:
        for i in groessen:
            dualvariablen[i][4].ub = 0.

    # Constraint (32b), left-hand-side expr
    inner_exprs = []
    for i in groessen:
        inner_expr = 1 * dualvariablen[i][0]\
                             - 1 * dualvariablen[i][1]\
                             - ret_time_plus[i] * dualvariablen[i][2]\
                             + ret_time_minus[i] * dualvariablen[i][3]\
                             - varianz_schranke[i] * dualvariablen[i][4]\
                             - zeit_diskret * factor * gp.quicksum(schlauch_rtd[i][j] * schlauchvariable[i][j] for j in range(anzahl_prozess))
        inner_exprs.append(inner_expr)

    reinheit_dual = gp.quicksum(inner_exprs[i]
                             for i in groessen)

    # params.fix requires purity as objective - in all other cases, the purity is bounded from below.
    # (32b). The right-hand side is 0 for the safe approximation and
    # -sum_s Delta_N^s for the upper bound of Corollary "enclosure".
    purity_constr = None
    if not params.fix:
        purity_constr = m.addConstr(reinheit_dual >= params.purity_rhs, "32b")

    constr_32c = []
    constr_32d = []
    for i in groessen:
        if i == params.wunschgroesse:
            fract_list = fract_p
        else:
            fract_list = fract_n
        zweimu = ret_time_plus[i] + ret_time_minus[i]

        zeile_c = []
        zeile_d = []
        for t in range(anzahl_prozess):
            # (32c)
            zeile_c.append(m.addConstr(1/factor * (a_s[i] * q0[i] *
                        fract_list[t]
                        - dualvariablen[i][0]
                        + dualvariablen[i][1]
                        + dualvariablen[i][2] * time_points[t]
                        - dualvariablen[i][3] * time_points[t]
                        + dualvariablen[i][4] * (time_points[t] ** 2 - zweimu * time_points[t]))
                        + schlauchvariable[i][t] >= 0, f'32c_{i}_{t}'))
        for t in range(anzahl_prozess - 1):
            # (32d)
            zeile_d.append(m.addConstr(1/factor * (a_s[i] * q0[i] *
                        fract_list[t]
                        - dualvariablen[i][0]
                        + dualvariablen[i][1]
                        + dualvariablen[i][2] * time_points[t+1]
                        - dualvariablen[i][3] * time_points[t+1]
                        + dualvariablen[i][4] * (time_points[t+1] ** 2 - zweimu * time_points[t+1]))
                        + schlauchvariable[i][t] >= 0, f'32d_{i}_{t}'))
        constr_32c.append(zeile_c)
        constr_32d.append(zeile_d)

    ##########################
    # OBJECTIVE ##############
    ##########################

    # lenght of fract interval
    zielfunktionMALTE = gp.quicksum(- i * jump_variable[i][0] + i * jump_variable[i][1] for i in range(anzahl_prozess + 1))
    # nominal fractionation volume
    zielfunktionMIT = (zeit_diskret / inst.totm_desired) * gp.quicksum(fractionierung[i] * matrix_nom[params.wunschgroesse][i] for i in range(anzahl_prozess + 1))
    # robust fractionation volume
    zielfunktionROBUST = 1 * dualvariablen[params.wunschgroesse][0]\
                             - 1 * dualvariablen[params.wunschgroesse][1]\
                             - ret_time_plus[params.wunschgroesse] * dualvariablen[params.wunschgroesse][2]\
                             + ret_time_minus[params.wunschgroesse] * dualvariablen[params.wunschgroesse][3]\
                             - varianz_schranke[params.wunschgroesse] * dualvariablen[params.wunschgroesse][4]\
                             - zeit_diskret * factor * gp.quicksum(schlauch_rtd[params.wunschgroesse][j] * schlauchvariable[params.wunschgroesse][j] for j in range(anzahl_prozess))

    # purity
    zielfunktionEVAL = sum(inner_expr for inner_expr in inner_exprs)

    if params.single_objective:
        # only the length of the fractionation interval. Its optimal value is the
        # same as the one of the objective hierarchy below, whose highest
        # priority is this objective.
        m.setObjective(zielfunktionMALTE, gp.GRB.MAXIMIZE)
        handles = {
            'fractionierung': fractionierung, 'fract_p': fract_p, 'fract_n': fract_n,
            'dualvariablen': dualvariablen, 'schlauchvariable': schlauchvariable,
            'jump_variable': jump_variable, 'inner_exprs': inner_exprs,
            'purity_constr': purity_constr, 'constr_32c': constr_32c, 'constr_32d': constr_32d,
            'obj_robust': zielfunktionROBUST, 'obj_nom': zielfunktionMIT,
            'obj_interval': zielfunktionMALTE, 'obj_eval': zielfunktionEVAL,
        }
        return m, handles

    # Objective expressions. Change priority for changing optimization goals.
    m.setObjectiveN(
        zielfunktionROBUST,
        index=0,
        priority=2,
        weight=-1.0,
        name="obj_robust"
    )

    m.setObjectiveN(
        zielfunktionMIT,
        index=1,
        priority=3,
        weight=-1.0,
        name="obj_nom"
    )
    m.setObjectiveN(
        zielfunktionMALTE,
        index=2,
        priority=4,
        weight=-1.0,
        name="obj_malte"
    )
    if params.fix:
        m.setObjectiveN(
            zielfunktionEVAL,
            index=3,
            priority=5,
            weight=-1.0,
            name="obj_eval"
        )

    handles = {
        'fractionierung': fractionierung,
        'fract_p': fract_p,
        'fract_n': fract_n,
        'dualvariablen': dualvariablen,
        'schlauchvariable': schlauchvariable,
        'jump_variable': jump_variable,
        'inner_exprs': inner_exprs,
        'purity_constr': purity_constr,
        'constr_32c': constr_32c,
        'constr_32d': constr_32d,
        'obj_robust': zielfunktionROBUST,
        'obj_nom': zielfunktionMIT,
        'obj_interval': zielfunktionMALTE,
        'obj_eval': zielfunktionEVAL,
    }
    return m, handles


def set_solver_parameters(m, params: Params):
    """Solver settings. Unchanged since version 0.1.1 except for the log file."""
    m.params.FeasibilityTol = 1e-9  # default is 1e-6
    #m.params.ScaleFlag = 1
    m.setParam("NumericFocus", 3)  # Highest level of numerical precision
    m.setParam("IntFeasTol", 1e-9)
    m.setParam('MIPGap', 0.00)
    if params.log_file is not None:
        m.setParam('LogFile', params.log_file)
    if params.time_limit is not None:
        m.setParam('TimeLimit', params.time_limit)


def solve_dro_model(time_points, matrix_nom_roh, matrix_min_roh, matrix_max_roh, params: Params):
    """Build and solve the MIP, print a summary and return the result.

    The returned dictionary contains, in addition to the entries of version
    0.1.1 (``optimal_frac``, ``running_time``), the solver statistics required
    for the model size table and the solution in terms of grid indices and
    times.
    """
    # start runtime tracking
    start = time.time()

    inst = baue_instanz(time_points, matrix_nom_roh, matrix_min_roh, matrix_max_roh, params)
    m, h = build_dro_model(inst, params)
    set_solver_parameters(m, params)

    # model statistics before solving (Gurobi counts presolve-independent sizes)
    m.update()
    stats = {
        'num_vars': m.NumVars,
        'num_bin_vars': m.NumBinVars,
        'num_int_vars': m.NumIntVars,
        'num_continuous_vars': m.NumVars - m.NumIntVars,
        'num_constrs': m.NumConstrs,
        'num_nz': m.NumNZs,
    }

    # optimize
    m.optimize()

    # uncomment to save lp/sol file
    #m.write("test.lp")
    #m.write("test.sol")

    # save end of runtime
    ende = time.time()
    status = m.status
    # a run stopped by the time limit can still have an incumbent
    feasible = m.SolCount >= 1
    if status == gp.GRB.Status.OPTIMAL:
        print("Feasible.")
    elif feasible:
        print(f"Feasible, but not solved to optimality (status {_status_name(status)}).")
    else:
        print("Infeasible.")

    stats.update({
        'status': int(status),
        'status_name': _status_name(status),
        'runtime': m.Runtime,
        'work': _attr(m, 'Work'),
        'node_count': m.NodeCount,
        # Gurobi does not expose MIPGap for multi-objective models; the runs are
        # terminated with MIPGap = 0, so an OPTIMAL status means gap 0.
        'mip_gap': _attr(m, 'MIPGap') if feasible else None,
        'obj_bound': _attr(m, 'ObjBound'),
        'sol_count': m.SolCount,
        'iter_count': m.IterCount,
        'wall_time': ende - start,
        'proven_optimal': status == gp.GRB.Status.OPTIMAL,
    })
    if params.single_objective and stats['obj_bound'] is not None:
        # the objective counts grid cells; convert the bound to minutes
        stats['length_bound'] = stats['obj_bound'] * inst.zeit_diskret

    result = {
        'stats': stats,
        'delta_N': inst.zeit_diskret,
        'params': {k: v for k, v in vars(params).items()},
        'running_time': ende - start,
    }

    if not feasible:
        result['optimal_frac'] = None
        print('Process length: ', inst.anzahl_prozess + 1)
        print('Runtime (seconds): ', ende - start)
        if params.stats_file:
            write_stats(result, params.stats_file)
        if params.keep_model:
            result['model'], result['handles'], result['instanz'] = m, h, inst
        return result

    # save fract values
    fraktionierung_werte = []
    for i in h['fractionierung']:
        fraktionierung_werte.append(int(round(i.x)))

    # analyze solution and determine begin and end of fractionation
    begin = end = None
    for i in range(inst.anzahl_prozess + 1):
        if h['jump_variable'][i][0].X >= 0.5:
            begin = i
        if h['jump_variable'][i][1].X >= 0.5:
            end = i

    inner_exprs = h['inner_exprs']
    a_s = inst.a_s
    purity_num = inner_exprs[params.wunschgroesse].getValue() / a_s[params.wunschgroesse]
    purity_den = sum(inner_exprs[i].getValue() / a_s[i] for i in inst.groessen)

    # some informative output
    print('Area:', aux.flaeche(time_points, matrix_nom_roh[params.wunschgroesse]))
    print('Process length: ', inst.anzahl_prozess + 1)
    print('Runtime (seconds): ', ende - start)
    print('OBJ nominal:', h['obj_nom'].getValue())
    print('OBJ robust:', h['obj_robust'].getValue())
    print('OBJ interval:', h['obj_interval'].getValue())
    if params.fix: print('OBJ eval:', h['obj_eval'].getValue())
    print(f"Fract interval indexes: ({begin},{end})")
    print(f"Fract interval timesteps: ({time_points[begin]},{time_points[end]})")
    if purity_den >= 1e-06:
        print(f"Worst Case Purity: {purity_num / purity_den}")
    else:
        print(f"Worst Case Purity: 0/0")
    print(f"Model size: {stats['num_continuous_vars']} continuous, "
          f"{stats['num_bin_vars']} binary, {stats['num_constrs']} constraints")
    print(f"Solver: runtime {stats['runtime']:.2f}s, {stats['node_count']:.0f} B&B nodes, "
          f"MIPGap {stats['mip_gap'] if stats['mip_gap'] is None else format(stats['mip_gap'], '.2e')}")

    result.update({
        'optimal_frac': fraktionierung_werte,
        'index_lower': begin,
        'index_upper': end,
        'x_lower': time_points[begin],
        'x_upper': time_points[end],
        'interval_length': time_points[end] - time_points[begin],
        'obj_nom': h['obj_nom'].getValue(),
        'obj_robust': h['obj_robust'].getValue(),
        'obj_interval': h['obj_interval'].getValue(),
        'purity_dual': (purity_num / purity_den) if purity_den >= 1e-06 else None,
        'inner_values': [inner_exprs[i].getValue() for i in inst.groessen],
        'a_coeff': [inst.a_coeff(i) for i in inst.groessen],
    })
    if params.fix:
        result['obj_eval'] = h['obj_eval'].getValue()

    if params.stats_file:
        write_stats(result, params.stats_file)

    if params.keep_model:
        result['model'], result['handles'], result['instanz'] = m, h, inst

    # return result
    return result


def _attr(m, name):
    """Model attribute or None if it is not available (e.g. MIPGap for
    multi-objective models)."""
    try:
        return getattr(m, name)
    except AttributeError:
        return None


def _status_name(status):
    for name in dir(gp.GRB.Status):
        if not name.startswith('_') and getattr(gp.GRB.Status, name) == status:
            return name
    return str(status)


def write_stats(result, path):
    serialisable = {k: v for k, v in result.items()
                    if k not in ('model', 'handles', 'instanz')}
    with open(path, 'w') as f:
        json.dump(serialisable, f, indent=2, default=float)
