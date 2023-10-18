# -*- coding: utf-8 -*-

import numpy as np


#%% Get surface tension of LJTS fluid with correlation by Vrabec et al., Molecular Physics 104 (2006), DOI: 10.1080/00268970600556774
def gamma_vrabec2006(T):
    '''
    Get surface tension of LJTS fluid (Vrabec 2006)

    :param float T: Temperature
    :return: float gamma: Surface tension
    '''
    
    # Equation 15
    T_c = 1.0779
    a = 2.08
    b = 1.21
    gamma = a*(1-(T/T_c))**b
    return gamma


#%% Get saturated densities and pressure by Vrabec et al., Molecular Physics 104 (2006), DOI: 10.1080/00268970600556774.
def sat_vrabec2006(T, radius=None):
    '''
    Get saturated densities by Vrabec et al., Molecular Physics 104 (2006). Equation numbers refer this paper.

    :param float T: Temperature
    :param float radius: reduced radius of spherical interface
    :return: float rhol, float rhov, float pv: Saturated liquid and vapor density and pressure
    '''
    
    tc,rc=1.0779,0.3190
    dt=tc-T
    a,b,c=0.5649,0.1314,0.0413
    rhol=rc+a*dt**(1/3.)+b*dt+c*dt**(3/2.)                      # equation 4
    a,b,c=0.5649,0.2128,0.0702
    rhov=rc-a*dt**(1/3.)+b*dt+c*dt**(3/2.)
    return rhol,rhov

if __name__ == "__main__":
    print('Running test with T = 0.8 ...')
    T = 0.8
    rhol_vrabec,rhov_vrabec,pv_vrabec = sat_vrabec2006(T)
    rhol_thol,rhov_thol,pv_thol = sat_thol2015(T)
    
    print('Sat. liquid dens. (Vrabec2006): '+str(rhol_vrabec))
    print('Sat. vapor dens. (Vrabec2006):  '+str(rhov_vrabec))
    print('Vapor pressure (Vrabec2006):  '+str(pv_vrabec))
    
    print('Sat. liquid dens. (Thol2015): '+str(rhol_thol))
    print('Sat. vapor dens. (Thol2015):  '+str(rhov_thol))
    print('Vapor pressure (Thol2015):  '+str(pv_thol))
    
    print('gamma Vrabec2006: '+str(gamma_vrabec2006(T)))
    
