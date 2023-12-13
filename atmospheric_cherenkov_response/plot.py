import numpy as np
import matplotlib
from matplotlib import colors
from matplotlib.colors import ListedColormap


def make_linear_cmap(color):
    o = np.ones(256)
    alpha = np.linspace(0, 1, 256)
    mycolors = np.c_[o * color[0], o * color[1], o * color[2], alpha]
    cmap = matplotlib.colors.ListedColormap(mycolors)
    return cmap


def particles():
    c = {}
    c["gamma"] = {}
    c["gamma"]["color"] = "black"
    c["gamma"]["rgb"] = (0, 0, 0)
    c["gamma"]["cmap"] = make_linear_cmap(c["gamma"]["rgb"])

    c["electron"] = {}
    c["electron"]["color"] = "blue"
    c["electron"]["rgb"] = (0, 0, 1)
    c["electron"]["cmap"] = make_linear_cmap(c["electron"]["rgb"])

    c["proton"] = {}
    c["proton"]["color"] = "red"
    c["proton"]["rgb"] = (1, 0, 0)
    c["proton"]["cmap"] = make_linear_cmap(c["proton"]["rgb"])

    c["helium"] = {}
    c["helium"]["color"] = "orange"
    c["helium"]["rgb"] = (1, 0.65, 0)
    c["helium"]["cmap"] = make_linear_cmap(c["helium"]["rgb"])
    return c


def energy_cmap(energy_start_GeV, energy_stop_GeV, cmap="nipy_spectral"):
    norm = matplotlib.colors.LogNorm(
        vmin=energy_start_GeV, vmax=energy_stop_GeV
    )
    return matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)


def config():
    c = {}
    c["particles"] = particles()
    c["label_unit_seperator"] = "$\\,/\\,$"
    c["rcParams"] = {"mathtext.fontset": "cm", "font.family": "STIXGeneral"}
    return c
