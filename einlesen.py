# general imports
import os

# parameter file
from params import Params

# generate folder file names
def rd4(params: Params):
    if params.sample == "small":
        daten = [
            'daten/Small_E',
            'daten/Small_E_min',
            'daten/Small_E_max',
        ]
        name="Small_E.txt"
    elif params.sample == "medium":
        daten = [
            'daten/Medium_E',
            'daten/Medium_E_min',
            'daten/Medium_E_max',
        ]
        name="Medium_E.txt"
    elif params.sample == "large":
        daten = [
            'daten/Large_E',
            'daten/Large_E_min',
            'daten/Large_E_max',
        ]
        name="Large_E.txt"
    else:
        assert False, "Invalid value for params.sample."

    return daten, name

# process the folder file names
def inputdata(ordner_nom, ordner_min, ordner_max):
    time_points, matrix_nom = _ein_ordner_einlesen(ordner_nom)
    _, matrix_min = _ein_ordner_einlesen(ordner_min)
    _, matrix_max = _ein_ordner_einlesen(ordner_max)
    return time_points, matrix_nom, matrix_min, matrix_max

# read one folder
def _ein_ordner_einlesen(ordnername):

    time_points = []
    peaks_matrix = []
    time_defined = False

    for path, subdirs, files in os.walk(ordnername):
        sor_files = sorted(files)
        for name in sor_files:
            peak = []
            datei = os.path.join(path, name)
            with open(datei, "r") as f:
                for line in f:
                    pair = line.split(' ')
                    if time_defined == False:
                        time_points.append(float(pair[0]))
                    peak.append(float(pair[1]))
            time_defined = True
            peaks_matrix.append(peak)
    return time_points, peaks_matrix
