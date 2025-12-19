# general imports
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

# auxiliary functions
import hilfsfunktionen as aux
# parameter file
from params import Params

def plot(data, optimal_frac, params: Params):
    # calculate envelopes
    env = aux.schlauch_chromatogramm(data[1], data[2], data[3])
    times = data[0]
    nominal = data[1]

    # transform
    nominal_peak = [max(profile) for profile in nominal]
    env_peak = [max(profile) for profile in env]

    # scale optimal_frac
    max_ref = max(max(nominal_peak), max(env_peak))
    max_opt = max(optimal_frac) if len(optimal_frac) > 0 else 0.0

    if max_opt > 0:
        scaled_optimal = [v * max_ref / max_opt for v in optimal_frac]
    else:
        scaled_optimal = optimal_frac

    fig, ax = plt.subplots()

    # plot nominal peaks
    for nomindex, nomi in enumerate(nominal):
        if nomindex == params.wunschgroesse:
            color = "red"
        else:
            color = "0.6"
        ax.plot(times, nomi, color=color, linestyle="-", label="")

    # plot envelope
    for envindex, enve in enumerate(env):
        if envindex == params.wunschgroesse:
            color = "firebrick"
        else:
            color = "black"

        ax.plot(times, enve, color=color, linestyle="-", label="")

    # optimal_frac area
    so = np.array(scaled_optimal)
    mask = so != 0
    idx = np.where(mask)[0]
    blocks = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
    for block in blocks:
        if len(block) == 0:
            continue
        x1 = times[block[0]]
        x2 = times[block[-1]]
        ax.axvspan(x1, x2, color="blue", alpha=0.3)

    # generate legend
    legend_elements = [
        Line2D([0], [0], color='black', lw=2, label='Envelope'),
        Line2D([0], [0], color='0.6', lw=2, label='Nominal'),
        Line2D([0], [0], color='blue', lw=4, alpha=0.3, label='Fractionation'),
        Line2D([0], [0], color='red', lw=2, label='Envelope Desired'),
        Line2D([0], [0], color='firebrick', lw=2, label='Nominal Desired')
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    # some axis settings and labels
    ax.set_xlabel("Time in min")
    ax.set_ylabel("Signal in pA")
    ax.set_ylim(-1, 16)
    fig.tight_layout()

    # save plot
    os.makedirs("output", exist_ok=True)
    fig.savefig("output/plot.pdf", dpi=300, bbox_inches="tight")

    # show plot
    plt.show()

    return fig, ax
