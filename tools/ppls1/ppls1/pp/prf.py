# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 13:14:26 2018

@author: mheinen
"""

import numpy as np

def remove_duplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def split_comp_prf(data):
    for di in data:
        header=di['header']
        tmpCount=len(header)-1
        header=remove_duplicates([s.strip('[0123456789]') for s in header])  # strip component id, remove duplicates
        di['header']=header
        numQ=len(header)-1
        di['numQ']=numQ
        numC=tmpCount//numQ-1
        di['numC']=numC
        
        # split data (quantities) componentwise 
        for qi in range(numQ):
            di[header[0]]=di['data'][:,0]  # extract coord vector
            di[header[1+qi]]=[ di['data'][:,(1+qi+ci*numQ)] for ci in range(numC+1) ]
            
def calc_Tcorr(dl_scal, dl_vect, writefreq):
    for d_scal,d_vect in zip(dl_scal,dl_vect):
        d_scal['Tcorr']=[]
        d_vect['Tycorr']=[]
        for ci in range(len(d_scal['Ekin_total'])):
            Ekin_total=np.multiply(d_scal['Ekin_total'][ci],writefreq)
#            T=np.multiply(Ekin_total,2)
            DOF_total=d_scal['DOF_total'][ci]
#            T=np.divide(T,DOF_total)
            N=np.divide(d_scal['DOF_trans'][ci],3)
            Ekin_drift=np.power(d_vect['vy'][ci],2)
            Ekin_drift=np.multiply(Ekin_drift,N)
            Ekin_drift=np.multiply(Ekin_drift,0.5)
            Ekin_T=np.subtract(Ekin_total,Ekin_drift)
            Tcorr=np.multiply(Ekin_T,2)
            Tcorr=np.divide(Tcorr,DOF_total)
            d_scal['Tcorr'].append(Tcorr)
            
            # Ty (Tz)
            Ekin_trans_y=np.multiply(d_vect['Ekin_trans,y'][ci],writefreq)
            Ty=d_vect['Ty'][ci]
            DOF_trans_y=np.divide(2*Ekin_trans_y,Ty)
            Ekin_Ty=np.subtract(Ekin_trans_y,Ekin_drift)
            Tycorr=np.multiply(Ekin_Ty,2)
            Tycorr=np.divide(Tycorr,DOF_trans_y)
            d_vect['Tycorr'].append(Tycorr)
    
        
    