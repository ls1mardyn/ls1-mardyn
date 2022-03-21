# -*- coding: utf-8 -*-

import numpy as np


#%% Get saturated densities by Vrabec et al., Molecular Physics 104 (2006).
def rho_vrabec2006(T):
    '''
    Get saturated densities by Vrabec et al., Molecular Physics 104 (2006)

    :param float T: Temperature
    :return: float rhol, float rhov: Saturated liquid and vapor density
    '''
    
    tc,rc=1.0779,0.3190
    dt=tc-t
    a,b,c=0.5649,0.1314,0.0413
    rhol=rc+a*dt**(1/3.)+b*dt+c*dt**(3/2.)
    a,b,c=0.5649,0.2128,0.0702
    rhov=rc-a*dt**(1/3.)+b*dt+c*dt**(3/2.)
    return rhol,rhov

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
    rhol,rhov = rho_vrabec2006(T)
    
    print('Sat. liquid dens. (Vrabec2006): '+str(rhol))
    print('Sat. vapor dens. (Vrabec2006):  '+str(rhov))

