import numpy as np # numpy
import petspy      # petsEOS python wrapper

# choose from a variety of input variables, e.g.
# DMASS (density) and T (temperature) or
# DMASS (density) and p (pressure) ...
# a comprehensive list can be found in pets.h, the nomenclature is taken from CoolProp (http://www.coolprop.org/)
# NB: Only for density and temperature as independent variables, the PeTS EOS can be evaluated directly. 
#     If density and any other quantity (pressure, internal energy, ...) is chosen, Newton algorithms are used to invert the EOS
integer_x = 12     # density
integer_y = 19     # temperature

# choose the values for the input variables
var_x = 0.5
var_y = 1.45

n_z = 7 # number of quantities to calculate
integer_z=np.zeros((n_z,1))

integer_z[0] = 15  # pressure
integer_z[1] = 21  # specific internal energy
integer_z[2] = 14  # specific enthalpy
integer_z[3] = 51  # specific Gibbs free energy
integer_z[4] = 18  # specific entropy
integer_z[5] = 32  # isochoric heat capacity
integer_z[6] = 34  # isobaric heat capacity

var_z=petspy.petseos(integer_x,var_x,integer_y,var_y,integer_z)


print(var_z[0])
print(var_z[1])
print(var_z[2])
print(var_z[3])
print(var_z[4])
print(var_z[5])
print(var_z[6])
