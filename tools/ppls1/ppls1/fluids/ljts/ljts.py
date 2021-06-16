# -*- coding: utf-8 -*-
"""
Created on Sun May 31 12:32:59 2020

@author: mheinen/homes
"""

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

if __name__ == "__main__":
    print('Running test with T = 0.8 ...')
    T = 0.8
    rhol,rhov = rho_vrabec2006(T)
    
    print('Sat. liquid dens. (Vrabec2006): '+str(rhol))
    print('Sat. vapor dens. (Vrabec2006):  '+str(rhov))

