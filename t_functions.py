"""
Function for Turbet atmospheric model

    - Declaration of constants required for the fits presented in the paper
    - Implementation of implicit functions
    - Implementation of the functions presented in the paper
"""

##########################################################################################

import numpy as np

##########################################################################################
# Model constants
                    # Values for surface temperature (eq. 1) #
k1 = 2.688
k2 = 1.099
k3 = 7.664e-1
k4 = 1.019
k5 = 4.683e-1
k6 = 4.224e-1

c1 = 3.401
c2 = 1.501e-1
c3 = -3.146e-2
c4 = 4.702e-2
c5 = -4.911e-3
c6 = 8.519e-3
c7 = -1.467e-2
c8 = -7.091e-3
c9 = -7.627e-3
c10 = 8.348e-3

                    # Values for eff temperature (eq. C.2) #
ek1 = 3.550
ek2 = 1.099
ek3 = 7.664e-1
ek4 = 1.310
ek5 = 4.683e-1
ek6 = 4.224e-1

ec1 = 2.846
ec2 = 1.555e-1
ec3 = 8.777e-2
ec4 = 6.045e-2
ec5 = 1.143e-2
ec6 = 1.736e-2
ec7 = 1.859e-2
ec8 = 4.314e-2
ec9 = 3.393e-2
ec10 = -1.034e-2

                    # General constants for conversions #

lim_pres_min = 2.7              # Minimum pressure validity     (bar)
lim_pres_max = 27e3             # Maximum pressure validity     (bar)
lim_S_min = 1                   # Minimum insolation validity     (earth unit)
lim_S_max = 30                  # Maximum insolation validity     (earth unit)

earth_m = 5.97216787e24              # Earth mass    (kg)    # CHANGE HERE WAS 5.8
earth_r = 6.3781e6              # Earth radius  (m)
earth_g = 9.81                  # Earth gravity (m.s-2)

conv_bar = 1e-5                 # Bar conversion (1 Pa = 1e-5 bar)
Gcst = 6.6743e-11                 # Gravitational constant (m3.kg-1.s-2)
Rcst = 8.31446262                    # Thermodynamical constant (J.K-1.mol-1)
wcp = 22064000                  # Water critcal pressure (Pa)
wct = 647.1                     # Water critical temperature (K)
wmm = 1.8e-2                    # Water molar mass (kg.mol-1)

P_transit = 2000        # Atmosphere cutoff (Pa)

##########################################################################################
"""
All inputs of functions have to be in earth units (M, R, g, S,...) and are converted within functions. All pressures have to be given in Pa, and are converted to bar within the funtions when necessary.
All outputs of functions are in earth units (and Pa).

Notation:
M - total mass
M_core - mass of core
M_atmo - mass of atmosphere

R - total radius
R_core - radius of core
Z_atmo - radius of atmosphere

S - Insolation of planet (in earth unit)
x_h2o - fraction of water (in mass compared to M, between 0 and 1)

P_surf - pressure of surface in Pa
T_surf - temperature at the surface of the core (K)
T_eff - temperature of the isothermal atmosphere (K)
"""
##########################################################################################
# Implicit functions
                    # Compute g from mass and radius #
def compute_g(M_i, R_i):
    M = earth_m*M_i
    R = earth_r*R_i
    return Gcst*M/R/R/earth_g

                    # Compute the mass of the atmosphere #
def compute_Matmo(M_core_i, x_h2o):
    M_core = earth_m*M_core_i
    return ((M_core*x_h2o)/(1-x_h2o))/earth_m

                    # Compute the pressure at the base of the atmosphere #
def compute_P_surf(M_atmo_i, g_i, R_i):
    M_atmo = earth_m*M_atmo_i
    g = earth_g*g_i
    R = earth_r*R_i
    return (M_atmo*g)/(4*np.pi*(R**2))      # Dimensional analysis gives Pa as output

##########################################################################################
# Model functions
"""
First set of equations is meant to compute the equation 1 in the paper, which leads to T surface. They use the first set of constants, conversion of units done in the main eq.

This set of equations use Pressure (P), gravity at the surface of the core planete (g) and the insolation of the planet (S).
"""
                    # Compute the X function of surface temperature #
def X(P):
    return (np.log10(P) - k1)/k2

                    # Compute the Y function of surface temperature #
def Y(g):
    return (np.log10(g) - k3)/k4

                    # Compute the Z function of surface temperature #
def Z(S):
    return (np.log10(S) - k5)/k6

                    # Computes the surface temperature (eq. 1) #
def T_surf(P_i, g_i, S):
    P = conv_bar*P_i
    g = earth_g*g_i
    x = X(P)
    y = Y(g)
    z = Z(S)
    return 10**( c1 + c2*x + c3*y + c4*z + c5*x*x + c6*x*y + c7*y*y + c8*z*z + c9*y*y*y + c10*z*z*z )

"""
Second set of equations is meant to compute the equation C.2 in the paper, which leads to T eff. They use the second set of constants, conversion of units done in the main eq.

This set of equations use water fraction (x_h2o), gravity at the surface of the core planete (g) and the insolation of the planet (S).
"""
                    # Compute the X function of eff temperature #
def eX(x_h2o):
    return (np.log10(x_h2o) - ek1)/ek2

                    # Compute the Y function of eff temperature #
def eY(g):
    return (g - ek3)/ek4

                    # Compute the Z function of eff temperature #
def eZ(S):
    return (np.log10(S) - ek5)/ek6

                    # Computes the eff temperature (eq. C.2) #
def T_eff(x_h2o, g_i, S):
    g = earth_g*g_i
    x = eX(x_h2o)
    y = eY(g)
    z = eZ(S)
    return 10**( ec1 + ec2*x + ec3*y + ec4*z + ec5*x*y + ec6*y*y + ec7*x*x*x + ec8*x*x*y + ec9*x*y*y + ec10*y*y*y*y )

"""
The last equation is dedicated to computing the atmosphere height, Z_atmo. The output is in earth unit once again.

T is the temperature (T_eff should be used, but might work with T_surf)
"""
                    # Computes the atmosphere height (eq. C.1) #
def Z_atmo(M_i, R_i, x_h2o, g_i, T):
    M = earth_m*M_i
    R = earth_r*R_i
    g = earth_g*g_i
    # Equation is large, so it's computed in multiple terms
    l1 = x_h2o/(1 - x_h2o)
    l2 = g*g/(4*np.pi*Gcst*P_transit)
    l3 = Rcst*T/(wmm*g)
    b1 = 1/(np.log10(l1*l2)*l3)
    b2 = 1/(b1 - (1/R))
    return b2/earth_r

##########################################################################################
# Other functions
"""
Semi-empirical law from Zeng et al. 2016, to compute R based on M and the core mass fraction. The radius computed corresponds to a planet without atmosphere.
The mass doesn't include the atmosphere mass; mass and radius are in earth units.
"""
                    # Empirical law from Zeng et al. 2016 #
def emp_R(M, cmf):
    return (1.07 - 0.21*cmf)*np.power(M, 1/3.7)
