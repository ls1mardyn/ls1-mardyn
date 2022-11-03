#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 23:48:38 2021

@author: mheinen/amartyads
"""
import sys
import os

CUR_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(CUR_DIR + '/../'))

from ppls1.ppls1.fluids.ljts.ljts import rho_vrabec2006

# For details, see Vrabec et al., Molecular Physics 104 (2006).
# The following variables are delared with respect to the default droplet experiment
# They may be edited as required

sig = 3.3916                #sigma value for Argon (see Table 3 in paper)
eps_times_kb = 137.9        #epsilon value for Argon multiplied with the Boltzmann constant (see Table 3 in paper)
scenario_temp = 110         #temperature of scenario, in K
drop_diameter = 1000        #droplet diameter, in Angstrom

Re=(drop_diameter/2.)/sig       #calculate reduced radius of droplet
t=scenario_temp/eps_times_kb    #calculate reduced temperature
rhol,rhov=rho_vrabec2006(t, Re) #get reduced densities
rhol=(rhol/(sig**3))            #convert rho into non reduced
rhov=(rhov/(sig**3))

# We now have densities in number of particles (or unified atomic mass) per Angstrom cubed
# This can be converted if need be, but is directly used in the given scenario
# For example, this can be converted to moles per litre by dividing with 6.02214e-4 (Avogadro's number, adjusted for converting Angstrom cubed to litre)

print('Droplet diameter: ' + str(drop_diameter) + '\nLiquid density: ' + str(rhol) + '\nVapour density: ' +str(rhov))

