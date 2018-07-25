# -*- coding: utf-8 -*-
"""
Created on Thu May 17 21:59:23 2018

@author: mheinen
"""

import numpy as np

#def imp_prf(fplist, flist=False):
def imp_prf(fplist):
    res=[]
#    fplist=[files]
#    if flist==True:
#        fplist=np.genfromtxt(files,dtype=np.str)
#        
    for fp in fplist:
        data=np.genfromtxt(fp,skip_header=1)
        numRows=data.shape[0]
        header=np.genfromtxt(fp,dtype=str,skip_footer=numRows)
        fpath=fp.rpartition('/')[0]
        fname=fp.rpartition('/')[-1]
        fext=fp.rpartition('.')[-1]
        finfo={'path':fpath,'name':fname,'ext':fext}
        ds={'header':header,
            'data':data,
            'finfo':finfo,
            'timestep':int(finfo['name'].rpartition('TS')[-1].split('.')[0])}
        res.append(ds.copy() )
    return res