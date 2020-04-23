"""
Grapher for results computed for the Turbet atmospheric model

    - Each graph is callable from 'main.py'
"""

##########################################################################################

import os
import numpy as np
from importlib import reload
import matplotlib.pyplot as pl

import t_functions as fct
reload(fct)

path_model = os.getcwd() + "/"
path_outputs = path_model + 'outputs/'

pl.style.use('seaborn-paper')

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

    pl.figure()

    print("    <i> Plotting core...")
    pl.plot(M, R, color='black', label='Core planet')
    print("    <i> Plotting atmospheres...")
    for x in Z:
        pl.plot(M, R + Z[x], label='$X_w$={0}%'.format(x))

    pl.legend(loc='best')
    pl.title("Mass radius relations with atmospheres ($P_t$={0} Pa)".format(fct.P_transit))
    pl.xlabel("Mass ($M_e$)")
    pl.ylabel("Radius ($R_e$)")

    #pl.ylim([0.0,3.0])
    pl.ylim([0.6, 1.3])     # Turbet limit on fig 2
    pl.xlim([0.15,2])

    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + 'atmospheres.png')
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

    pl.figure()

    print("    <i> Plotting relation...")
    pl.plot(M, g, color='black')

    pl.title("Gravity of core at the surface")
    pl.ylabel("Gravity (g)")
    pl.xlabel("Mass ($M_e$)")
    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + 'gravity.png')
    pl.close()
    print("> Finished graphing gravity!")
    print("--------------------------------------------------------------")

def surface_pressure(M, P):
    print("> Graphing surface pressure computed for all the planets...")
    print("    <i> Applying cutoffs: {0} - {1} Me...".format(cutoff[0], cutoff[1]))
    md = np.where(M<cutoff[0])
    Md = np.where(M>cutoff[1])
    M = np.delete(M, md)
    M = np.delete(M, Md)
    for i in P:
        P[i] = np.delete(P[i], md)
        P[i] = np.delete(P[i], Md)

    pl.figure()

    print("    <i> Plotting relations...")
    for x in P:
        pl.plot(M, P[x], label='$X_w$={0}%'.format(x))

    pl.title("Pressure at the surface")
    pl.xlabel("Mass ($M_e$)")
    pl.ylabel("Pressure (Pa)")
    pl.legend(loc='best')
    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + "surface_pressure.png")
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

    pl.figure()

    print("    <i> Plotting relations...")
    for x in T:
        pl.plot(g, T[x], label='$X_w$={0}%'.format(x))

    pl.title("Surface temperature")
    pl.xlabel("Surface gravity (g)")
    pl.ylabel("Temperature (K)")
    pl.legend(loc='best')
    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + "temperature_surface.png")
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

    print("    <i> Plotting relations...")
    for x in T:
        pl.plot(g, T[x], label='$X_w$={0}%'.format(x))

    pl.title("Effective temperature")
    pl.xlabel("Surface gravity (g)")
    pl.ylabel("Temperature (K)")
    pl.legend(loc='best')

    print("    <i> Saving to outputs folder...")
    pl.savefig(path_outputs + "temperature_effective.png")
    pl.close()
    print("> Finished graphing effective temperature!")
    print("--------------------------------------------------------------")
