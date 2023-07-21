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
    rhov=rc-a*dt**(1/3.)+b*dt+c*dt**(3/2.)                      # equation 5
    a,b,c=3.1664,5.9809,0.01498
    pv = np.exp(a - b/T + c/(T**(3)))
    
    if radius != None:
        a,b,c,d,e=0.1485,0.0471,1.272,0.4944,5.0090
        rhol=rhol+(a-b*dt**2)/(radius+c-d*dt**-1+e*dt)  # equation 18
        
        a,b,c,d,e=0.0541,0.1682,1.0810,0.2490,10.87
        rhov=rhov+(-a+b*T**2)/(radius-c-d*dt**-1+e*dt)      # equation 19
        
    return rhol,rhov,pv


#%% Get saturated densities by Thol et al., Int J Thermophys 36 (2015), DOI: 10.1007/s10765-014-1764-4.
def sat_thol2015(T):
    '''
    Get saturated densities by Thol et al., Int J Thermophys 36 (2015)

    :param float T: Temperature
    :return: float rhol, float rhov: Saturated liquid and vapor density
    '''
    
    coeff = {"pv" : np.array([-6.21,1.5,-1.92,2.2,-4.76]),
             "rhol": np.array([1.45,-0.172,-0.298,0.295]),
             "rhov": np.array([1.59809,-0.09975,-0.4774,-2.33736])}
    exp   = {"pv" : np.array([1,1.5,3.25,4.85,6.63]),
             "rhol" : np.array([0.334,0.667,1.25,1.92]),
             "rhov" : np.array([1,1.5,5.94,0.41452])}
    
    crit1  = {"pv" : 1.086 , "rhol" : 1.086 , "rhov" : 1.086}
    crit2  = {"pv" : 0.101 , "rhol" : 0.319 , "rhov" : 0.319}

    summe = sum(coeff['rhol']*(1.0 - T/crit1['rhol'])**exp['rhol'])
    rhol = crit2['rhol']*np.exp(summe)
    
    summe = sum(coeff['rhov']*(1.0 - T/crit1['rhov'])**exp['rhov'])
    rhov = crit2['rhov']*(1.0 + summe)
    
    summe = sum(coeff['pv']*(1.0 - T/crit1['pv'])**exp['pv'])
    pv = crit2['pv']*np.exp(crit1['pv']/T*summe)
    
    return rhol,rhov,pv


#%% Convert Gibbs free energy from PeTS to ms2 value
def g_PeTS2ms2(g_pets, T, rho):
    '''
    Convert Gibbs free energy from PeTS to ms2 value
    
    :param float g_pets: Gibbs free energy from PeTS
    :param float T: Temperature
    :param float rho: Density
    :return: float g_ms2: Gibbs free energy from ms2
    '''
    
    rhocrit = 0.319
    Tcrit   = 1.086
    delta0 = (0.001/0.8)/rhocrit
    tau0   = Tcrit/0.8
    ig1 = -2.5/tau0
    ig2 = 1.5 - np.log(delta0) - 1.5*np.log(tau0)
    tau   = Tcrit/T
    delta = rho/rhocrit
    alphaId = np.log(delta) + 1.5*np.log(tau) + ig1*tau + ig2
    g_res = g_pets/T - 1.0 - alphaId
    g_ms2 = g_res + np.log(rho)
    return g_ms2


#%% Convert Gibbs free energy from ms2 to PeTS value
def g_ms22PeTS(g_ms2, T, rho):
    '''
    Convert Gibbs free energy from ms2 to PeTS value
    
    :param float g_ms2: Gibbs free energy from ms2
    :param float T: Temperature
    :param float rho: Density
    :return: float g_pets: Gibbs free energy from PeTS
    '''
    
    rhocrit = 0.319
    Tcrit   = 1.086
    delta0 = (0.001/0.8)/rhocrit
    tau0   = Tcrit/0.8
    ig1 = -2.5/tau0
    ig2 = 1.5 - np.log(delta0) - 1.5*np.log(tau0)
    tau   = Tcrit/T
    delta = rho/rhocrit
    alphaId = np.log(delta) + 1.5*np.log(tau) + ig1*tau + ig2
    g_res = g_ms2 - np.log(rho)
    g_pets = T*(1.0 + alphaId + g_res)
    return g_pets


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
    
