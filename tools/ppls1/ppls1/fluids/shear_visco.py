# -*- coding: utf-8 -*-
"""
Created on 2022

@author: homes
"""

import numpy as np
import pandas as pd

#%% Get dynamic viscosity with correlation by Lemmon/Refprop (International Journal of Thermophysics, Vol. 25, No. 1, January 2004)
def eta_lemmon(T,rho,fluid,units='reduced'):
    '''
    Get dynamic viscosity (Lemmon2004)

    :param float T: Temperature (SI: in [K])
    :param float rho: Density (SI: in [mol/l])
    :param str fluid: Name of fluid (Air|Argon|Nitrogen|Oxygen|LJTS|LJfull)
    :param str units: System of input/output units (reduced|SI; default: reduced)
    :return: float eta: dynamic viscosity (SI: in [kg/(m s)])
    '''
    
    def Omega(Tstar,fluid):
        summe=0
        for i in range(0,len(tab2)):
            bi=tab2['bi'][i]
            summe+=(bi*np.log(Tstar)**i)
        return np.exp(summe)

    def eta0(T,fluid):
        M,eps_kb,sig=tab1[fluid][3],tab1[fluid][4],tab1[fluid][5]
        #print("M,eps_kb,sig=",M,eps_kb,sig)
        Tstar=T/eps_kb
        eta0=(0.0266958*np.sqrt(M*T)) / (sig**2*Omega(Tstar,fluid))
        #eta0=(0.168729283*np.sqrt(T)) / (sig**2*Omega(Tstar,fluid))  # Monika
        return eta0 #8.18940 #22.7241 #8.18940 #eta0

    def etar(tau,delta,fluid):
        summe=0
        for ri in range(0,len(tab3[fluid])):
            i,Ni,ti,di,li=tab3[fluid]['i'][ri],tab3[fluid]['Ni'][ri],tab3[fluid]['ti'][ri],tab3[fluid]['di'][ri],tab3[fluid]['li'][ri]
            if li==0:
                gamma=0
            else:
                gamma=1
            summe+=Ni*tau**ti*delta**di*np.exp(-gamma*delta**li)
        return summe

    def eta(T,rho,fluid):
        Tc,rhoc=tab1[fluid][0],tab1[fluid][1]
        tau,delta=Tc/T,rho/rhoc
        eta=eta0(T,fluid)+etar(tau,delta,fluid)
        return eta


    # tab1
    d={'Parameter': ['T_c (K)','r_c (mol*dm^−3)','p_c (MPa)','M (g*mol^−1 )','eps/kb (K)','sig (nm)','psi_0 (nm)','Tau','q_D (nm)','T_ref (K)'],
    'Nitrogen': [126.192,11.1839,3.3958,28.01348,98.94,0.3656,0.17,0.055,0.40,252.384],
    'Argon': [150.687,13.40743,4.863,39.948,143.2,0.335,0.13,0.055,0.32,301.374],
    'Oxygen': [154.581,13.63,5.043,31.9988,118.5,0.3428,0.24,0.055,0.51,309.162],
    'Air': [132.6312,10.4477,3.78502,28.9586,103.3,0.360,0.11,0.055,0.31,265.262]}
    tab1=pd.DataFrame(data=d)

    # tab2
    d={'i': [0,1,2,3,4],
    'bi': [0.431,-0.4623,0.08406,0.005341,-0.00331]}
    tab2=pd.DataFrame(data=d)

    # tab3
    tab3={'Nitrogen':None,'Argon':None,'Oxygen':None,'Air':None}
    d={'i': [1,2,3,4,5],
    'Ni': [10.72,0.03989,0.001208,-7.402,4.620],
    'ti': [0.1,0.25,3.2,0.9,0.3],
    'di': [2,10,12,2,1],
    'li': [0,1,1,2,3]}
    tab3['Nitrogen']=pd.DataFrame(data=d)

    d={'i': [1,2,3,4,5,6],
    'Ni': [12.19,13.99,0.005027,-18.93,-6.698,-3.827],
    'ti': [0.42,0.0,0.95,0.5,0.9,0.8],
    'di': [1,2,10,5,1,2],
    'li': [0,0,0,2,4,4]}
    tab3['Argon']=pd.DataFrame(data=d)

    d={'i': [1,2,3,4,5],
    'Ni': [17.67,0.4042,0.0001077,0.3510,-13.67],
    'ti': [0.05,0.0,2.10,0.0,0.5],
    'di': [1,5,12,8,1],
    'li': [0,0,0,1,2]}
    tab3['Oxygen']=pd.DataFrame(data=d)

    d={'i': [1,2,3,4,5],
    'Ni': [10.72,1.122,0.002019,-8.876,-0.02916],
    'ti': [0.2,0.05,2.4,0.6,3.6],
    'di': [1,4,9,1,8],
    'li': [0,0,0,1,1]}
    tab3['Air']=pd.DataFrame(data=d)

    # tab4
    tab4={'Nitrogen':None,'Argon':None,'Oxygen':None,'Air':None}
    d={'i': [1,2,3,4,5,6,7,8,9],
    'Ni': [1.511,2.117,-3.332,8.862,31.11,-73.13,20.03,-0.7096,0.2672],
    'ti': [None,-1.0,-0.7,0.0,0.03,0.2,0.8,0.6,1.9],
    'di': [None,None,None,1,2,3,4,8,10],
    'li': [None,None,None,0,0,1,2,2,2]}
    tab4['Nitrogen']=pd.DataFrame(data=d)

    d={'i': [1,2,3,4,5,6,7,8,9,10],
    'Ni': [0.8158,-0.4320,0.0,13.73,10.07,0.7375,-33.96,20.47,-2.274,-3.973],
    'ti': [None,-0.77,-1.0,0.0,0.0,0.0,0.8,1.2,0.8,0.5],
    'di': [None,None,None,1,2,4,5,6,9,1],
    'li': [None,None,None,0,0,0,2,2,2,4]}
    tab4['Argon']=pd.DataFrame(data=d)

    d={'i': [1,2,3,4,5,6,7,8,9],
    'Ni': [1.036,6.283,-4.262,15.31,8.898,-0.7336,6.728,-4.374,-0.4747],
    'ti': [None,-0.9,-0.6,0.0,0.0,0.3,4.3,0.5,1.8],
    'di': [None,None,None,1,3,4,5,7,10],
    'li': [None,None,None,0,0,0,2,2,2]}
    tab4['Oxygen']=pd.DataFrame(data=d)

    d={'i': [1,2,3,4,5,6,7,8,9],
    'Ni': [1.308,1.405,-1.036,8.743,14.76,-16.62,3.793,-6.142,-0.3778],
    'ti': [None,-1.1,-0.3,0.1,0.0,0.5,2.7,0.3,1.3],
    'di': [None,None,None,1,2,3,7,7,11],
    'li': [None,None,None,0,0,2,2,2,2]}
    tab4['Air']=pd.DataFrame(data=d)
    
    allFluids = ['Nitrogen','Argon','Oxygen','Air','LJTS','LJfull']
    
    if fluid not in allFluids:
        print('Fluid not yet supported')
        print('Possible fluids:')
        print(allFluids)
        return 0.0
    
    na=6.02214076e23
    kb=1.380649e-23
    u_mass=1.660539e-27
    
    if fluid == 'LJTS':
        fluid = 'Argon'
        units = 'reduced'
        sig=3.3916
        eps=137.9
        mass=39.948
        tc=1.08
        # rc=0.31
    elif fluid == 'LJfull':
        fluid = 'Argon'
        units = 'reduced'
        sig=3.3952
        eps=116.79
        mass=39.948
        tc=1.32
        # rc=0.31
    else:
        units = 'SI'
        sig=1.0
        eps=1.0
        mass=1.0
        tc=150.687
    
    if units == 'reduced':
        T = T/tc*150.687
        rho = rho/sig**3*1e30/na*1e-3
        # refTime, refLambda: tref, lref
        tref=sig*1e-10*np.sqrt(mass*u_mass/(kb*eps))
        eta_ref=1e6*mass*u_mass/(sig*1e-10*tref) # Lemmon gives [uPa s] (see Table V in Paper)
    elif units == 'SI':
        pass
    else:
        print('Unit unknown')
        return 0.0
    
    eta_value=eta(T, rho, fluid)
    if units == 'reduced':
        return eta_value/eta_ref
    elif units == 'SI':
        return eta_value


#%% Get dynamic viscosity of LJTS fluid with correlation by Lautenschläger2019
def eta_lauten(T,rho):
    '''
    Get dynamic viscosity of LJTS fluid (Lautenschläger)

    :param float T: Temperature
    :param float rho: Density
    :return: float eta_l: dynamic viscosity
    '''
    
    b = 0.5736
    c = 1.38
    d = 2.438
    coeff = np.matrix('0.1786,  0, 0; 0.06529, 1, 0; -1.158,  0, 1; 1.347, 0, 2; 0.3194, -1, 1')
    eta_l = b*rho*np.exp(c*rho+((d*(rho**4)-rho)/T))
    
    for i in range(5):
            eta_l = eta_l + coeff[i,0]*(T**(coeff[i,1]))*(rho**(coeff[i,2]))
    return eta_l


#%% Get dynamic viscosity of LJfull fluid with correlation by Galliero, Ind. Eng. Chem. Res., vol. 44, 2005
def eta_galliero(T,rho):
    '''
    Get dynamic viscosity of LJfull fluid (Galliero)

    :param float T: Temperature
    :param float rho: Density
    :return: float eta_l: dynamic viscosity
    '''
    
    def omega_22(T):
        a_22=1.16145
        b_22=0.14874
        c_22=0.52487
        d_22=0.7732
        e_22=2.16178
        f_22=2.43787
        omega_22 = a_22/T**b_22+c_22*np.exp(-d_22*T)+e_22*np.exp(-f_22*T)
        return omega_22
    
    b1=0.062692
    b2=4.095577
    b3=-8.743269e-6
    b4=11.12492
    b5=2.542477e-6
    b6=14.863984

    tmp1 = 0.17630924*np.sqrt(T)/omega_22(T)
    tmp2 = b1*(np.exp(b2*rho)-1)+b3*(np.exp(b4*rho)-1)+b5*(np.exp(b6*rho)-1)/T**2

    return tmp1+tmp2


#%% Tests
if __name__ == '__main__':
    fluid = 'LJTS'
    T = 0.95
    rho = 0.7
    print(f'Running test with {fluid} ...')
    print('Eta Lemmon:         '+str(eta_lemmon(T,rho,fluid)))
    print('Eta Lautenschlager: '+str(eta_lauten(T,rho)))

    fluid = 'LJfull'
    T = 1.1
    rho = 0.7
    print(f'Running test with {fluid} ...')
    print('Eta Lemmon:   '+str(eta_lemmon(T,rho,fluid)))
    print('Eta Galliero: '+str(eta_galliero(T,rho)))

    fluid = 'Argon'
    T = 200.0
    rho = 10.0
    print(f'Running test with {fluid} ...')
    print('Eta Lemmon:     '+str(eta_lemmon(T,rho,fluid)))
    print('Eta Lemmon Lit: 25.5662')  # Value from Table V in Lemmon2004
