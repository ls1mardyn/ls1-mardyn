# -*- coding: utf-8 -*-
"""
Created on 2022

@author: homes
"""

import numpy as np

#%% Get shear viscosity of LJTS fluid with correlation by Lautenschläger2019
def eta_lauten(T,rho):
    '''
    Get shear viscosity of LJTS fluid (Lautenschläger)

    :param float T: Temperature
    :param float rho: Density
    :return: float eta_l: Shear viscosity
    '''
    
    b = 0.5736
    c = 1.38
    d = 2.438
    coeff = np.matrix('0.1786,  0, 0; 0.06529, 1, 0; -1.158,  0, 1; 1.347, 0, 2; 0.3194, -1, 1')
    eta_l = b*rho*np.exp(c*rho+((d*(rho**4)-rho)/T))
    for i in range(5):
            eta_l = eta_l + coeff[i,0]*(T**(coeff[i,1]))*(rho**(coeff[i,2]))
    return eta_l

if __name__ == '__main__':
    print('Running test with LJTS ...')
    fluid = 'LJTS'
    T = 0.8
    rho = 0.55

    print('Eta Lautenschlager: '+str(eta_lauten(T,rho)))

