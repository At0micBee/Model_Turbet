"""
Main file for Turbet atmospheric model

Computations done from 't_functions.py'.
Graphs produced from 't_grapher.py'.

Implementation by H.G. Vivien - hugo.vivien@pm.me
"""
##########################################################################################
# Loading required libraries
import os
import numpy as np
from importlib import reload

import t_functions as fct
import t_grapher as graph
reload(fct)
reload(graph)

##########################################################################################
# Defining path to work folders

path_model = os.getcwd() + "/"
path_src = path_model + 'sources/'
path_input = path_model + 'inputs/'
path_outputs = path_model + 'outputs/'

##########################################################################################
# Printing info about model properties
print("##############################################################")
print("           >>> Turbet model implementation (V.2) <<<          ")
print("##############################################################")
print("")
print("> Figures and data files produced are stored in:")
print("> {0}".format(path_outputs))
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("")

##########################################################################################
# Selecting data to be used
"""
The final result will be influenced by what type of composition is choosen for the core planet. This relation can be changed by plugging in a different file for MR relations, here we include the 'MR_rock.txt' file, which contains the MR relations computed by Zeng et al. 2016 for a pure rocky composition (MgSiO3).

!!! Format necessary: 2 columns, first is M, second is R, (no header) in earth units. !!!
"""
print("> General informations about current options:")
file = "MR_rock.txt"
planet = np.loadtxt(path_src + file)
print("    <i> File used for Mass-Radius relation: {0}".format(file))
print("    <i> M from {0} to {1} earth masses".format(planet[:,0][0], planet[:,0][-1]))

"""
The height of the atmosphere will slightly vary depending on the selected transition pressure to determine the limit of the atmosphere with space.
The default value (as stated in Turbet et al. 2019) is 0.1 Pa, but can be modified in the 't_functions.py' module if required (stored as 'P_transit').
"""
print("    <i> Atmosphere pressure cutoff selected: {0} Pa".format(fct.P_transit))

"""
The insolation of the planet is required to compute the temperatures, and can be set to values between 0 and 30 to respect the limit of the model. In some cases an array is used to study the variations, but the main value that is going to be used throughout the program is set here. Value is earth unit (1S = 1366 W.m-2).
"""
S = 1
print("    <i> Main solar insulation is set to: {0} earth insolation".format(S))

print("--------------------------------------------------------------")

##########################################################################################
# Basic computation of certain values
                    # Compute gravity associated to file loaded #
gravity = fct.compute_g(planet[:,0], planet[:,1])

                    # Prepare the water fraction properties #
"""
Water fraction values to be used to compute pressures in the model, enter values in percent. The computation part then uses those to get T_surf and T_eff.
"""
frac_water = np.array([0.01, 0.1, 1, 3, 5])
mass_atmosphere, pressure_water, temperature_surf, temperature_eff, Z_atm = {}, {}, {}, {}, {}
for val in frac_water:
    mass_atmosphere[val] = fct.compute_Matmo(planet[:,0], val/100)
    pressure_water[val] = fct.compute_P_surf(mass_atmosphere[val], gravity, planet[:,1])
    """
    for i in range(len(pressure_water[val])):
        pressure_water[val][i] = fct.wcp
    """
    temperature_surf[val] = fct.T_surf(pressure_water[val], gravity, S)
    temperature_eff[val] = fct.T_eff(val/100, gravity, S)
    Z_atm[val] = fct.Z_atmo(planet[:,0], planet[:,1], val/100, gravity, temperature_eff[val])

##########################################################################################
# Saving computed data
                    # Single file with combination of all data computed #
all = open(path_outputs + "results.dat", "w+")
all.write("wf\t\t\tM\t\t\tR\t\t\tg\t\t\tP\t\t\t\t\t\t\t\t\tT_surf\t\t\t\tT_eff\t\t\tZ_atmo\n")
for frac in frac_water:
    for i in range(len(planet[:,0])):
        all.write("{0:8.6f}\t{1:8.6f}\t{2:8.6f}\t{3:8.6f}\t{4:28.20f}\t{5:16.6f}\t{6:12.6f}\t{7:12.6f}\n".format(frac, planet[:,0][i], planet[:,1][i], gravity[i], pressure_water[frac][i], temperature_surf[val][i], temperature_eff[val][i], Z_atm[val][i]))
all.close()

##########################################################################################
# User inputs
                    # Calling the graphing functions #
graph.atmospheres(planet[:,0], planet[:,1], Z_atm)
graph.gravity_planet(planet[:,0], gravity)
graph.surface_pressure(planet[:,0], gravity, pressure_water)
graph.temp_surf(planet[:,0], gravity, temperature_surf)
graph.temp_eff(planet[:,0], gravity, temperature_eff)
