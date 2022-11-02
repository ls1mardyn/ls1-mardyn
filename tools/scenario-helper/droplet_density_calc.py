#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 23:48:38 2021

@author: mheinen
"""

# For details, see Vrabec et al., Molecular Physics 104 (2006).
# Equation numbers refer above paper
def rho_vrabec2006(t,Re=None):
    tc,rc=1.0779,0.3190
    dt=tc-t
    a,b,c=0.5649,0.1314,0.0413
    rhol=rc+a*dt**(1/3.)+b*dt+c*dt**(3/2.)                      # equation 4
    a,b,c=0.5649,0.2128,0.0702
    rhov=rc-a*dt**(1/3.)+b*dt+c*dt**(3/2.)                      # equation 5
    
    if Re != None:
        a,b,c,d,e=0.1485,0.0471,1.272,0.4944,5.0090
        rhol=rhol+(a-b*(tc-t)**2)/(Re+c-d*(tc-t)**-1+e*(tc-t))  # equation 18
        
        a,b,c,d,e=0.0541,0.1682,1.0810,0.2490,10.87
        rhov=rhov+(-a+b*t**2)/(Re-c-d*(tc-t)**-1+e*(tc-t))      # equation 19
    
    return rhol,rhov

# The following variables are delared with respect to the default droplet experiment
# They may be edited as required

sig = 3.3916                #sigma value for Argon
eps_times_kb = 137.9        #epsilon value for Argon multiplied with the Boltzmann constant
scenario_temp = 110         #temperature of scenario, in K
drop_diameter = 1000        #droplet diameter, in Angstrom

Re=(drop_diameter/2.)/sig       #calculate reduced radius of droplet
t=scenario_temp/eps_times_kb    #calculate reduced temperature
rhol,rhov=rho_vrabec2006(t,Re)
rhol=(rhol/(sig**3))            #convert rho into non reduced
rhov=(rhov/(sig**3))

# We now have densities in number of particles (or unified atomic mass) per Angstrom cubed
# This can be converted if need be, but is drectly used in the given scenario
# For example, this can be converted to moles per litre by dividing with 6.02214e-4 (Avogadro's number, adjusted for converting Angstrom cubed to litre)

print('Droplet diameter: ' + str(drop_diameter) + '\nLiquid density: ' + str(rhol) + '\nVapour density: ' +str(rhov))

