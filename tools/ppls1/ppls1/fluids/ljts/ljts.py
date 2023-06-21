# -*- coding: utf-8 -*-
"""
Created on Sun May 31 12:32:59 2020

@author: mheinen/homes
"""

#%% Get saturated densities by Vrabec et al., Molecular Physics 104 (2006).
def rho_vrabec2006(t,Re=None):
    '''
    Get saturated densities by Vrabec et al., Molecular Physics 104 (2006). Equation numbers refer this paper.

    :param float T: Temperature, float Re: reduced radius of spherical interface
    :return: float rhol, float rhov: Saturated liquid and vapor density
    '''
    
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

if __name__ == "__main__":
    print('Running test with T = 0.8 ...')
    T = 0.8
    rhol,rhov = rho_vrabec2006(T)
    
    print('Sat. liquid dens. (Vrabec2006): '+str(rhol))
    print('Sat. vapor dens. (Vrabec2006):  '+str(rhov))

