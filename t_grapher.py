"""
Grapher for results computed for the Turbet atmospheric model

    - Each graph is callable from 'main.py'
"""

##########################################################################################

import os
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as pl
from importlib import reload

import t_functions as fct
reload(fct)

path_model = os.getcwd() + "/"
path_outputs = path_model + 'outputs/'

pl.style.use('seaborn-paper')
pl.rc("axes", labelsize=12.45)
pl.rc("xtick", labelsize=9.5)
pl.rc("ytick", labelsize=9.5)
pl.rc("legend", fontsize=10)

##########################################################################################

cutoff = [0, 22]            # Mass limitation for plots

def atmospheres(M, R, Z):
    print("> Graphing mass-radius relations, and atmospheres...")
    print("    <i> Applying cutoffs: {0} - {1} Me...".format(cutoff[0], cutoff[1]))
    md = np.where(M<cutoff[0])
    Md = np.where(M>cutoff[1])
    M = np.delete(M, md)
    M = np.delete(M, Md)
    R = np.delete(R, md)
    R = np.delete(R, Md)
    for i in Z:
        Z[i] = np.delete(Z[i], md)
        Z[i] = np.delete(Z[i], Md)
        Z[i][np.where(Z[i]<0)] = np.nan

    print("    <i> Plotting core...")

    pl.figure()

    pl.plot(M, R, color='black', label='Core planet')

    pl.ylim([0.6, 1.3])     # Turbet limit on fig 2
    pl.xlim([0.0,2])

    pl.xlabel("Mass ($M_e$)")
    pl.ylabel("Radius ($R_e$)")

    print("    <i> Plotting atmospheres...")
    pl.twinx()
    pl.twiny()

    pl.plot(M, R, color='black', label='Core planet')       # Adding here too for legend

    for x in Z:
        pl.plot(M, R + Z[x], label='$X_w$={0}%'.format(x))

    pl.legend(loc='best')
    #pl.title("Mass radius relations with atmospheres ($P_t$={0} Pa)".format(fct.P_transit))

    #pl.ylim([0.0,3.0])
    pl.ylim([0.6, 1.3])     # Turbet limit on fig 2
    pl.xlim([0.0,2])

    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + 'atmospheres.pdf')
    pl.close()
    print("> Finished graphing mass-radius relations, and atmospheres!")
    print("--------------------------------------------------------------")

def gravity_planet(M, g):
    print("> Graphing gravity computed for all the planets...")
    print("    <i> Applying cutoffs: {0} - {1} Me...".format(cutoff[0], cutoff[1]))
    md = np.where(M<cutoff[0])
    Md = np.where(M>cutoff[1])
    M = np.delete(M, md)
    M = np.delete(M, Md)
    g = np.delete(g, md)
    g = np.delete(g, Md)

    pl.figure(figsize=(6.6, 4.95))
    pl.xscale("log")
    pl.yscale("log")

    print("    <i> Plotting relation...")
    pl.plot(M, g, color='black', label="Surface gravity")

    #pl.title("Gravity of core at the surface")
    pl.ylabel("Gravity ($g_e$)")
    pl.xlabel("Mass ($M_e$)")

    pl.legend(loc="best")

    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + 'gravity.pdf')
    pl.close()
    print("> Finished graphing gravity!")
    print("--------------------------------------------------------------")

def surface_pressure(M, g, P):
    print("> Graphing surface pressure computed for all the planets...")
    print("    <i> Applying cutoffs: {0} - {1} Me...".format(cutoff[0], cutoff[1]))
    md = np.where(M<cutoff[0])
    Md = np.where(M>cutoff[1])
    M = np.delete(M, md)
    M = np.delete(M, Md)
    g = np.delete(g, md)
    g = np.delete(g, Md)
    for i in P:
        P[i] = np.delete(P[i], md)
        P[i] = np.delete(P[i], Md)

    pl.figure(figsize=(6.6, 4.95))
    pl.xscale("log")
    pl.yscale("log")

    # Plotting ranges
    pl.fill_between(fct.lim_g, fct.lim_pres[0]/fct.conv_bar, fct.lim_pres[1]/fct.conv_bar, color="grey", alpha=0.5)

    print("    <i> Plotting relations...")
    for x in P:
        pl.plot(g, P[x], label='$X_w$={0}%'.format(x))

    #pl.title("Pressure at the surface")
    pl.xlabel("Gravity ($g_e$)")
    pl.ylabel("Pressure (Pa)")
    pl.legend(loc='best')

    pl.twinx()
    pl.yscale("log")

    for x in P:
        pl.plot(g, P[x]*fct.conv_bar, alpha=0.0)

    pl.ylabel("Pressure (bar)")

    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + "surface_pressure.pdf")
    pl.close()
    print("> Finished graphing surface pressure!")
    print("--------------------------------------------------------------")

def temp_surf(M, g, T):
    print("> Graphing surface temperature computed for all planets...")
    print("    <i> Applying cutoffs: {0} - {1} Me...".format(cutoff[0], cutoff[1]))
    md = np.where(M<cutoff[0])
    Md = np.where(M>cutoff[1])
    M = np.delete(M, md)
    M = np.delete(M, Md)
    g = np.delete(g, md)
    g = np.delete(g, Md)
    for i in T:
        T[i] = np.delete(T[i], md)
        T[i] = np.delete(T[i], Md)

    pl.figure(figsize=(6.6, 4.95))
    pl.xscale("log")

    print("    <i> Plotting relations...")
    for x in T:
        pl.plot(g, T[x], label='$X_w$={0}%'.format(x))

    #pl.title("Surface temperature")
    pl.xlabel("Surface gravity ($g_e$)")
    pl.ylabel("Temperature (K)")
    pl.legend(loc='best')

    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + "temperature_surface.pdf")
    pl.close()
    print("> Finished graphing surface temperature!")
    print("--------------------------------------------------------------")

def temp_eff(M, g, T):
    print("> Graphing effective temperature computed for all planets...")
    print("    <i> Applying cutoffs: {0} - {1} Me...".format(cutoff[0], cutoff[1]))
    md = np.where(M<cutoff[0])
    Md = np.where(M>cutoff[1])
    M = np.delete(M, md)
    M = np.delete(M, Md)
    g = np.delete(g, md)
    g = np.delete(g, Md)
    for i in T:
        T[i] = np.delete(T[i], md)
        T[i] = np.delete(T[i], Md)

    pl.figure()

    # Plotting ranges
    pl.fill_between(fct.lim_g, fct.lim_T[0], fct.lim_T[1], color="grey", alpha=0.5)

    print("    <i> Plotting relations...")
    for x in T:
        pl.plot(g, T[x], label='$X_w$={0}%'.format(x))

    #pl.title("Effective temperature")
    pl.xlabel("Surface gravity (g)")
    pl.ylabel("Temperature (K)")
    pl.legend(loc='best')

    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + "temperature_effective.pdf")
    pl.close()
    print("> Finished graphing effective temperature!")
    print("--------------------------------------------------------------")
