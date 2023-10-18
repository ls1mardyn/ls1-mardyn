# -*- coding: utf-8 -*-
"""
Created on 2020

@author: homes
"""

import numpy as np
import pandas as pd

#%% Get thermal conductivity with correlation by Lemmon/Refprop (International Journal of Thermophysics, Vol. 25, No. 1, January 2004)
def lambda_lemmon(T,rho,fluid,units='reduced'):
    '''
    Get thermal conductivity (Lemmon2004)

    :param float T: Temperature (SI: in [K])
    :param float rho: Density (SI: in [mol/l])
    :param str fluid: Name of fluid (Air|Argon|Nitrogen|Oxygen|LJTS|LJfull)
    :param str units: System of input/output units (reduced|SI; default: reduced)
    :return: float lam: Thermal conductivity (SI: in [mW/(m K)])
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

    def lam0(T,fluid):
        Tc=tab1[fluid][0]
        tau=Tc/T
        N1,N2,N3=tab4[fluid]['Ni'][0],tab4[fluid]['Ni'][1],tab4[fluid]['Ni'][2]
        t2,t3=tab4[fluid]['ti'][1],tab4[fluid]['ti'][2]
        term1=N1*eta0(T,fluid)
        term2=N2*tau**t2
        term3=N3*tau**t3
        return term1+term2+term3

    def lamr(tau,delta,fluid):
        summe=0
        for ri in range(3,len(tab4[fluid])):
            i,Ni,ti,di,li=tab4[fluid]['i'][ri],tab4[fluid]['Ni'][ri],tab4[fluid]['ti'][ri],tab4[fluid]['di'][ri],tab4[fluid]['li'][ri]
            if li==0:
                gamma=0
            else:
                gamma=1
            summe+=Ni*tau**ti*delta**di*np.exp(-gamma*delta**li)
        return summe

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
        tref=sig*1e-10*np.sqrt(mass*1.660538e-27/(kb*eps))
        lref=1e3*kb/(sig*1e-10*tref) # [mW/m-K]
    elif units == 'SI':
        pass
    else:
        print('Unit unknown')
        return 0.0
    
    Tc,rhoc=tab1[fluid][0],tab1[fluid][1]
    tau,delta=Tc/T,rho/rhoc
    lam=lam0(T,fluid)+lamr(tau,delta,fluid)
    if units == 'reduced':
        return lam/lref
    elif units == 'SI':
        return lam

#%% Get thermal conductivity of LJTS fluid with correlation by Lautenschläger2018
def lambda_lauten(T,rho):
    '''
    Get thermal conductivity of LJTS fluid (Lautenschläger)

    :param float T: Temperature
    :param float rho: Density
    :return: float lambda_l: Thermal conductivity
    '''
    
    coeff = np.matrix('-0.08633,4.377,-14.25,25.49,-5.527; 0.5723,0.6753,0.7847,-0.03917,0.;-0.127,-0.1697,-0.04694,0.,0.;0.02022,0.01147,0.,0.,0.;-0.001088,0.,0.,0.,0.')
    lambda_l = 0
    #print(coeff)
    for i in range(5):
        for j in range((4-(i-1))):
            lambda_l = lambda_l + coeff[i,j]*(T**(i))*(rho**(j))
    return lambda_l

#%% Get thermal conductivity of 2CLJQ fluid with data by Dissertation Fernandez (DATA MIGHT BE WRONG!)
def lambda_fernandez(T,rho,L,QQ):
    '''
    Get thermal conductivity of 2CLJQ fluids (Fernandez) (DATA MIGHT BE WRONG)

    :param float T: Temperature
    :param float rho: Density
    :param float L: Elongation
    :param float QQ: Quadrupole momentum
    :return: float lambda_f: Thermal conductivity
    '''
    
    coeffAll = pd.DataFrame(index=(0,1,2,3,4),columns=(0.0,0.2,0.4,0.505,0.6,0.8))
    coeffAll[0.0][0] = np.matrix('-350.3841, 856.1682, -639.2193, 162.5303; 106.6538, -185.0944, 75.2985, 0; -9.9704, 9.1858, 0, 0; 0.29199, 0, 0, 0')
    coeffAll[0.2][0] = np.matrix('-21.9622, -873.4945, 1908.4555, -1039.6442; 146.9692, -135.3694, 9.0697, 0; -26.2684, 14.9843, 0, 0; 1.4343, 0, 0, 0')
    coeffAll[0.4][0] = np.matrix('-587.4527, 1893.1875, -2234.3074, 979.3953; 293.6705, -559.1689, 303.1685, 0; -53.9842, 44.7712, 0, 0; 3.6539, 0, 0, 0')
    coeffAll[0.505][0] = np.matrix('-875.2489, 2075.498, -908.8419, -382.0348; 656.4768, -1266.9356, 554.6385, 0; -145.138, 147.7971, 0, 0; 10.3225, 0, 0, 0')
    coeffAll[0.6][0] = np.matrix('-78.5074, 1111.4156, -2770.4706, 2034.3935; -102.441, 64.5953, 101.8584, 0; 39.5456, -34.3996, 0, 0; -3.7003, 0, 0, 0')
    coeffAll[0.8][0] = np.matrix('714.2697, -3772.0523, 6678.6368, -3947.2282; -440.5246, 1514.2072, -1266.3172, 0; 95.0351, -166.3094, 0, 0; -6.7309, 0, 0, 0')
    coeffAll[0.0][1] = np.matrix('-1478.5665, 3907.3411, -4206.1729, 1760.6173; 338.3333, -363.5331, 90.9423, 0; -42.2544, 23.3625, 0, 0; 1.7318, 0, 0, 0')
    coeffAll[0.2][1] = np.matrix('294.2808, -1751.1416, 2786.4317, -1400.7331; 32.542, 49.8367, -47.1321, 0; -10.0898, -0.62067, 0, 0; 0.78447, 0, 0, 0')
    coeffAll[0.4][1] = np.matrix('-534.3498, 2504.7658, -4177.6214, 2437.3733; 138.196, -316.9759, 211.5714, 0; -22.1081, 20.3718, 0, 0; 1.428, 0, 0, 0')
    coeffAll[0.505][1] = np.matrix('439.6406, -1455.6178, 1246.3454, 25; -263.6754, 717.0805, -488.4614, 0; 40.9023, -55.2715, 0, 0; -2.1202, 0, 0, 0')
    coeffAll[0.6][1] = np.matrix('-639.2652, 3158.3732, -5610.2553, 3622.3223; 266.2187, -707.0956, 481.7547, 0; -52.1631, 66.5667, 0, 0; 3.5846, 0, 0, 0')
    coeffAll[0.8][1] = np.matrix('-1059.0601, 5487.4754, -9189.2428, 5002.2977; 645.7857, -2307.0005, 2033.4437, 0; -123.7106, 222.8679, 0, 0; 7.8035, 0, 0, 0')
    coeffAll[0.0][2] = np.matrix('-1361.1972, 3003.0695, -2288.9532, 602.4054; 379.1246, -532.7006, 202.112, 0; -36.3821, 23.5023, 0, 0; 1.2462, 0, 0, 0')
    coeffAll[0.2][2] = np.matrix('-1301.9582, 4342.6756, -4898.8609, 1869.6092; 306.9072, -629.8384, 348.3126, 0; -28.8755, 24.6601, 0, 0; 1.1664, 0, 0, 0')
    coeffAll[0.4][2] = np.matrix('869.4924, -3866.3119, 5742.8487, -2868.5558; -200.2602, 589.4388, -395.5782, 0; 15.2965, -28.9655, 0, 0; -0.013995, 0, 0, 0')
    coeffAll[0.505][2] = np.matrix('1350.371, -4909.2547, 6591.7455, -3307.1483; -669.9114, 1389.3443, -703.1474, 0; 130.7081, -137.5478, 0, 0; -8.3662, 0, 0, 0')
    coeffAll[0.6][2] = np.matrix('-1366.1329, 5441.2814, -7063.798, 3049.3076; 729.9065, -1986.9296, 1324.507, 0; -124.5596, 174.0688, 0, 0; 6.8002, 0, 0, 0')
    coeffAll[0.8][2] = np.matrix('718.833, -2331.0111, 2150.3384, -223.3088; -663.3877, 1612.4201, -993.5434, 0; 183.7877, -221.299, 0, 0; -16.9466, 0, 0, 0')
    coeffAll[0.0][3] = np.matrix('-820.7873, 305.6523, 1142.3128, -692.5913; 403.7404, -405.2148, 72.3344, 0; -46.3092, 27.0691, 0, 0; 1.6025, 0, 0, 0')
    coeffAll[0.2][3] = np.matrix('-4486.5901, 12704.4308, -12911.9042, 4724.9086; 1189.1622, -1930.266, 822.9726, 0; -128.1281, 97.11, 0, 0; 4.9103, 0, 0, 0')
    coeffAll[0.4][3] = np.matrix('1239.5249, -4118.5454, 5295.154, -2472.7007; -476.5067, 805.9502, -387.199, 0; 81.2743, -60.4462, 0, 0; -5.0478, 0, 0, 0')
    coeffAll[0.505][3] = np.matrix('783.9174, -3902.7458, 6480.6104, -3654.1332; -217.8145, 679.3288, -460.211, 0; 27.0441, -51.8968, 0, 0; -0.65518, 0, 0, 0')
    coeffAll[0.6][3] = np.matrix('-715.5637, 2938.1326, -3785.5952, 1551.1147; 362.2669, -1052.0136, 742.9552, 0; -56.4999, 83.6003, 0, 0; 2.9318, 0, 0, 0')
    coeffAll[0.8][3] = np.matrix('-978.9972, 2801.9092, -2484.1852, 892.2952; 873.3364, -1693.0772, 656.0376, 0; -258.5682, 280.6847, 0, 0; 23.5794, 0, 0, 0')
    coeffAll[0.0][4] = np.matrix('2951.4429, -5494.95, 3452.5335, -740.1541; -771.3518, 949.1368, -285.3397, 0; 67.4283, -42.2227, 0, 0; -1.926, 0, 0, 0')
    coeffAll[0.2][4] = np.matrix('-870.0858, 2054.0956, -1552.4803, 318.9541; 254.4181, -401.0782, 191.9749, 0; -24.6808, 14.8274, 0, 0; 0.97277, 0, 0, 0')
    coeffAll[0.4][4] = np.matrix('801.3577, -3019.5806, 4078.0475, -1974.5161; -226.4268, 481.9491, -234.241, 0; 27.9996, -32.6628, 0, 0; -1.0393, 0, 0, 0')
    coeffAll[0.505][4] = np.matrix('1475.1695, -4502.7318, 5127.8273, -2153.2332; -774.1717, 1406.3979, -672.5891, 0; 147.6992, -128.8999, 0, 0; -9.5853, 0, 0, 0')
    coeffAll[0.6][4] = np.matrix('-895.2727, 1507.7097, 1118.6626, -2319.0262; 762.7905, -1513.2158, 794.1916, 0; -160.3815, 150.4297, 0, 0; 11.6694, 0, 0, 0')
    coeffAll[0.8][4] = np.matrix('-1343.4126, 6679.1253, -11396.6054, 6739.9362; 640.6699, -1960.1382, 1545.676, 0; -116.8778, 168.1148, 0, 0; 7.7285, 0, 0, 0')
    
    coeff = coeffAll[L][QQ]
    
    lambda_f = 0
    #print(coeff)
    for i in range(4):
        for j in range(4-i):
            lambda_f = lambda_f + coeff[i,j]*(T**(i))*(rho**(j))
    return lambda_f

#%% Get thermal conductivity of LJTS fluid with data/fit by Guevara/Homes
def lambda_guevara_homes(T,rho):
    '''
    Get thermal conductivity of LJTS fluid (Guevara/Homes)

    :param float T: Temperature
    :param float rho: Density
    :return: float lambda_l: Thermal conductivity
    '''
    
    p00=-364.1
    p10=-1003
    p01=2292
    p20=-21.11
    p11=2238
    p02=-3751
    p30=601.2
    p21=-1482
    p12=13.63
    p03=1494
    lambda_gh = p00 + p10*T + p01*rho + p20*T**2 + p11*T*rho + p02*rho**2 + p30*T**3 + p21*T**2*rho + p12*T*rho**2 + p03*rho**3
    return lambda_gh

if __name__ == '__main__':
    print('Running test with LJTS ...')
    fluid = 'LJTS'
    T = 0.8
    rho = 0.55

    print('Lambda Lemmon:         '+str(lambda_lemmon(T,rho,fluid)))
    print('Lambda Lautenschlager: '+str(lambda_lauten(T,rho)))
    print('Lambda Guevara/Homes:  '+str(lambda_guevara_homes(T,rho)))

