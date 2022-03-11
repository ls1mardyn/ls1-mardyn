import petspy      # petsEOS python wrapper
import argparse


# Copy or link libpets.so to /usr/local/lib/


# choose from a variety of input variables, e.g.
# DMASS (density) and T (temperature) or
# DMASS (density) and p (pressure) ...
# a comprehensive list can be found in pets.h, the nomenclature is taken from CoolProp (http://www.coolprop.org/)
# NB: Only for density and temperature as independent variables, the PeTS EOS can be evaluated directly. 
#     If density and any other quantity (pressure, internal energy, ...) is chosen, Newton algorithms are used to invert the EOS

variables = {'density':12,
            'enthalpy':14,
            'pressure':15,
            'entropy':18,
            'temperature':19,
            'internalEnergy':21,
            'isochoricHeatCapacity':32,
            'isobaricHeatCapacity':34,
            'gibbsEnergy':51,
            'chemicalPot':512}

arguments = {'temperature':None, 'density':None, 'pressure':None}

input_int = list()
input_val = list()

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbosity', action='count', default=0, help='Verbosity')
parser.add_argument('-T', '--temperature', type=float, help='Input: Temperature')
parser.add_argument('-d', '--density', type=float, help='Input: Density')
parser.add_argument('-p', '--pressure', type=float, help='Input: Pressure')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='Output, e.g.: temperature, density, pressure, internalEnergy, enthalpy, entropy, isochoricHeatCapacity, isobaricHeatCapacity, gibbsEnergy, chemicalPot \n'+
                    'Example: pets_args -T 1.22 -d 0.4 -o pressure')
args = parser.parse_args()
arguments['temperature'] = args.temperature
arguments['density']     = args.density
arguments['pressure']    = args.pressure
output = args.output

for key in arguments:
    if arguments[key] != None:
        input_int.append(variables[key])
        input_val.append(arguments[key])

if (len(input_int) > 2) or (len(input_val) > 2): print('WARNING: Too many arguments given!')

integer_z = variables[output]

var_z = petspy.petseos(input_int[0],input_val[0],input_int[1],input_val[1],integer_z)

if args.verbosity >= 1:
    print(f'{output} : {var_z[0]}')
else:
    print(var_z[0])


