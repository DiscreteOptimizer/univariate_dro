# general imports
import numpy as np
import time

# gurobi
import gurobipy as gp

# auxiliary functions
import hilfsfunktionen as aux
# params
from params import Params

def solve_dro_model(time_points, matrix_nom_roh, matrix_min_roh, matrix_max_roh, params: Params):
    # start runtime tracking
    start = time.time()

    # scaling factors to increase numeric stability
    factor = 1e06
    q_factor = 100.

    # initialize particle mass vector
    q0 = np.zeros(len(matrix_nom_roh))
    # calculate mass of particle i
    for i in range(len(matrix_nom_roh)):
        q0[i] = aux.flaeche(time_points, matrix_nom_roh[i])

    # initialize scaled data
    matrix_nom = np.zeros((len(matrix_nom_roh),len(matrix_nom_roh[0])))
    matrix_min = np.zeros((len(matrix_nom_roh),len(matrix_nom_roh[0])))
    matrix_max = np.zeros((len(matrix_nom_roh),len(matrix_nom_roh[0])))

    # scale data...
    for i, (el1, el2, el3, el4) in enumerate(zip(q0, matrix_nom_roh, matrix_min_roh, matrix_max_roh)):
        matrix_nom[i] = el2 / el1
        matrix_min[i] = el3 / el1
        matrix_max[i] = el4 / el1

    # particle masses
    for i in range(len(q0)):
        q0[i] *= q_factor

    # particle indices
    groessen = list(i for i in range(len(matrix_nom)))
    # number of time steps
    anzahl_prozess = len(time_points) - 1

    # total time interval length
    aux_sum = 0
    for i in range(len(time_points) - 1):
        aux_sum += time_points[i + 1] - time_points[i]

    # calculate \delta_N
    zeit_diskret = aux_sum / (len(time_points) - 1)  # berechnen des Zeitschritts als Mittelwert aller Zeitschritte

    # total mass of desired peak ist 1...
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

    # variance bound
    varianz_schranke =  [schwank_var_global * nomi - muplus * muminus + zeit_diskret ** 2 / 4. for nomi, muplus, muminus in zip(dict_var_nom, ret_time_plus, ret_time_minus)]

    # calculate envelopes
    schlauch_rtd = aux.schlauch_chromatogramm(matrix_nom, matrix_min, matrix_max)

    ##########################
    # MODEL ##################
    ##########################

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
    # (32b)
    if not params.fix: m.addConstr(reinheit_dual >= 0, "32b")

    a_s = []
    for i in groessen:
        if i == params.wunschgroesse:
            a_s.append(1 - params.reinheit)
        else:
            a_s.append(-params.reinheit)

    for i in groessen:
        if i == params.wunschgroesse:
            fract_list = fract_p
        else:
            fract_list = fract_n
        zweimu = ret_time_plus[i] + ret_time_minus[i]

        for t in range(anzahl_prozess):
            # (32c)
            m.addConstr(1/factor * (a_s[i] * q0[i] *
                        fract_list[t]
                        - dualvariablen[i][0]
                        + dualvariablen[i][1]
                        + dualvariablen[i][2] * time_points[t]
                        - dualvariablen[i][3] * time_points[t]
                        + dualvariablen[i][4] * (time_points[t] ** 2 - zweimu * time_points[t]))
                        + schlauchvariable[i][t] >= 0, f'32c_{i}_{t}')
        for t in range(anzahl_prozess - 1):
            # (32d)
            m.addConstr(1/factor * (a_s[i] * q0[i] *
                        fract_list[t]
                        - dualvariablen[i][0]
                        + dualvariablen[i][1]
                        + dualvariablen[i][2] * time_points[t+1]
                        - dualvariablen[i][3] * time_points[t+1]
                        + dualvariablen[i][4] * (time_points[t+1] ** 2 - zweimu * time_points[t+1]))
                        + schlauchvariable[i][t] >= 0, f'32d_{i}_{t}')

    ##########################
    # OBJECTIVE ##############
    ##########################

    # lenght of fract interval
    zielfunktionMALTE = gp.quicksum(- i * jump_variable[i][0] + i * jump_variable[i][1] for i in range(anzahl_prozess + 1))
    # nominal fractionation volume
    zielfunktionMIT = (zeit_diskret / totm_desired) * gp.quicksum(fractionierung[i] * matrix_nom[params.wunschgroesse][i] for i in range(anzahl_prozess + 1))
    # robust fractionation volume
    zielfunktionROBUST = 1 * dualvariablen[params.wunschgroesse][0]\
                             - 1 * dualvariablen[params.wunschgroesse][1]\
                             - ret_time_plus[params.wunschgroesse] * dualvariablen[params.wunschgroesse][2]\
                             + ret_time_minus[params.wunschgroesse] * dualvariablen[params.wunschgroesse][3]\
                             - varianz_schranke[params.wunschgroesse] * dualvariablen[params.wunschgroesse][4]\
                             - zeit_diskret * factor * gp.quicksum(schlauch_rtd[params.wunschgroesse][j] * schlauchvariable[params.wunschgroesse][j] for j in range(anzahl_prozess))

    # purity
    zielfunktionEVAL = sum(inner_expr for inner_expr in inner_exprs)

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


    # uncomment for advanced settings
    #m.params.FeasibilityTol = 1e-9 #default is 1e-6
    #m.params.ScaleFlag = 1
    #m.setParam("NumericFocus", 3)  # Highest level of numerical precision
    #m.setParam("IntFeasTol", 1e-9)
    #m.setParam('MIPGap', 0.05)

    # optimize
    m.optimize()

    # uncomment to save lp/sol file
    #m.write("test.lp")
    #m.write("test.sol")

    # save end of runtime
    ende = time.time()
    if m.status == gp.GRB.Status.OPTIMAL:
        print("Feasible.")
    else:
        print("Infeasible.")

    # save fract values
    fraktionierung_werte = []
    for i in fractionierung:
        fraktionierung_werte.append(int(i.x))

    # analyze solution and determine begin and end of fractionation
    for i in range(anzahl_prozess + 1):
        if jump_variable[i][0].X >= 0.5:
            begin = i
        if jump_variable[i][1].X >= 0.5:
            end = i

    # some informative output
    print('Area:', aux.flaeche(time_points, matrix_nom_roh[params.wunschgroesse]))
    print('Process length: ',len(fractionierung))
    print('Runtime (seconds): ', ende - start)
    print('OBJ nominal:', zielfunktionMIT.getValue())
    print('OBJ robust:', zielfunktionROBUST.getValue())
    print('OBJ interval:', zielfunktionMALTE.getValue())
    if params.fix: print('OBJ eval:', zielfunktionEVAL.getValue())
    print(f"Fract interval indexes: ({begin},{end})")
    print(f"Fract interval timesteps: ({time_points[begin]},{time_points[end]})")
    if sum(inner_exprs[i].getValue()/a_s[i] for i in groessen) >= 1e-06:
        print(f"Purity: {inner_exprs[params.wunschgroesse].getValue()/a_s[params.wunschgroesse]/sum(inner_exprs[i].getValue()/a_s[i] for i in groessen)}")
    else:
        print(f"Purity: 0/0")

    # store result
    result = {
        'optimal_frac': fraktionierung_werte,
        'running_time': ende - start
    }

    # return result
    return result
