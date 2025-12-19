# general imports
import numpy as np

# parameter file
from params import Params

# calculate expectations
def baue_mu_list(time_points, matrix):
    mu_list = []     # dictionary für Erwartungswerte erstellen
    for row in matrix:   # ueber alle Peaks (Listen) laufen
        mu_list.append(sum(element*time for time, element in zip(time_points, row)) / sum(row))
    return mu_list

# calculate variances
def baue_var_list(time_points, matrix, mu_list):
    var_list = []
    for row, mu in zip(matrix, mu_list):
        var_list.append(sum(time ** 2 * element for time, element in zip(time_points, row)) / sum(row) - mu ** 2)
    return var_list

# calculate particle mass
def flaeche(time_points, peak):
    flae = 0
    for i in range(len(time_points)-1):
        flae += (time_points[i+1] - time_points[i]) * peak[i]
    return flae

# re-calculate yield purity
def calculate_yieldpurity(data_short, fraktionierungsvektor, params: Params):
    time_points = data_short[0]
    matrix_nom = np.array(data_short[1])
    matrix_min = np.array(data_short[2])
    matrix_max = np.array(data_short[3])

    num_particles = len(matrix_nom)

    # indices where fract vector has ones
    idfrac = [i for i, val in enumerate(fraktionierungsvektor) if val == 1]

    # corresponding time points
    frac = [time_points[i] for i in idfrac]
    frac_mass = np.zeros((3, num_particles))

    for comp in range(num_particles):
        frac_mass[0, comp] = flaeche(frac, matrix_nom[comp][idfrac])
        frac_mass[1, comp] = flaeche(frac, matrix_min[comp][idfrac])
        frac_mass[2, comp] = flaeche(frac, matrix_max[comp][idfrac])

    total_mass = 0
    for comp in range(num_particles):
        if comp < params.wunschgroesse:
            total_mass += frac_mass[2, comp]
        elif comp > params.wunschgroesse:
            total_mass += frac_mass[1, comp]
        else:
            total_mass += frac_mass[0, comp]
    purity = frac_mass[0, params.wunschgroesse] / total_mass
    return purity

# aggregate density matrices
def aggregate_matrix(matrix, aggregation_factor, params: Params):
    new_matrix = []
    for index, row in enumerate(matrix): # partikelgroessen
        new_matrix.append([])
        part_sum = 0
        for index2, entry in enumerate(row):
            part_sum += entry
            if index2 % aggregation_factor == aggregation_factor-1:
                if part_sum / aggregation_factor >= 1e-05:
                    new_matrix[-1].append(part_sum/aggregation_factor)# + int(index == params.wunschgroesse) * 1e-05)
                else:
                    new_matrix[-1].append(0.)# + int(index == params.wunschgroesse) * 1e-05)
                part_sum = 0
    return new_matrix

# aggratate time step vectors
def aggregate_matrix_index(matrix, aggregation_factor):
    new_matrix = matrix[aggregation_factor-1::aggregation_factor]
    return new_matrix

find_max_index = lambda liste: max((wert, i) for i, wert in enumerate(liste))

# Envelope calculation
def schlauch_chromatogramm(ma1, ma2, ma3):
    matrix = [[0. for _ in row] for row in ma1]
    for row1, row2, row3, row in zip(ma1, ma2, ma3, matrix):  # Ueber alle Peaks gehen
        max1, ind1 = find_max_index(row1)     # Max Wert und Index mit Hilfsfunktion 'find_max_index' bestimmen
        max2, ind2 = find_max_index(row2)
        max3, ind3 = find_max_index(row3)
        for i, (el1, el2, el3) in enumerate(zip(row1, row2, row3)):    # ueber alle Zeitpunkte gehen
            if i < min(ind1, ind2, ind3) or i > max(ind1, ind2, ind3): # bestimmen ob, man zwischen den beiden aeußeren Peaks ist
                row[i] = max(el1, el2, el3)     # Wenn nicht, die Huelle um die Peaks "legen"
            else:
                row[i] = max(max1, max2, max3)  # Wenn ja, Werte auf maximal Wert setzten
    return matrix

# remove time steps where fract mass is 0 for all particles
def remove_zeros(data):
    my_list = [0] * len(data[0])
    for dat in data[1:]:
        for partsize in dat:
            for index, element in enumerate(partsize):
                my_list[index] += element
    first_index_nonzero = 0
    last_index_nonzero = 0
    for index, element in enumerate(my_list):
        if abs(element) >= 1e-09 and last_index_nonzero == 0:
            first_index_nonzero = index
        if abs(element) >= 1e-09:
            last_index_nonzero = index
    data[0] = data[0][first_index_nonzero:last_index_nonzero+1]
    for index1, dat in enumerate(data):
        if index1 ==  0: continue
        for index2, partsize in enumerate(dat):
            data[index1][index2] = data[index1][index2][first_index_nonzero:last_index_nonzero+1]
    return data
