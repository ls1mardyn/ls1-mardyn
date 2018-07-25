# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 11:29:05 2018

@author: mheinen
"""
import numpy as np
from struct import *

def imp_coal_bin(fplist, flist=False):
    res=[]
#    fplist=[files]
#    if flist==True:
#        fplist=np.genfromtxt(files,dtype=np.str,)
        
#    print(fplist)
#    print("len:",len(fplist))
        
    for fp in fplist:
        f=open(fp,'rb')
        try:
            numBins=(unpack('<L',f.read(4)))[0]
            numShells=(unpack('<L',f.read(4)))[0]
            dt=np.dtype('<f8')
            data=np.fromfile(f,dtype=dt)
            data=data.reshape(numShells,numBins)
            fpath=fp.rpartition('/')[0]
            fname=fp.rpartition('/')[-1]
            fext=fp.rpartition('.')[-1]
            finfo={'path':fpath,'name':fname,'ext':fext}
            ds={'numBins':numBins,
                'numShells':numShells,
                'data':data,
                'finfo':finfo,
                'timestep':int(finfo['name'].rpartition('TS')[-1].split('.')[0])}
            res.append(ds.copy() )
        finally:
            f.close()
    return res

def imp_coal_txt(files, flist=False):
    res=[]
    fplist=[files]
    if flist==True:
        fplist=np.genfromtxt(files,dtype=np.str,)
        
#    print(fplist)
#    print("len:",len(fplist))
        
    for fp in fplist:
        dt=np.dtype('<f8')
        data=np.genfromtxt(fp,dt)
        numShells,numBins=data.shape
        fpath=fp.rpartition('/')[0]
        fname=fp.rpartition('/')[-1]
        fext=fp.rpartition('.')[-1]
        finfo={'path':fpath,'name':fname,'ext':fext}
        ds={'numBins':numBins,
            'numShells':numShells,
            'data':data,
            'finfo':finfo,
            'timestep':int(finfo['name'].rpartition('TS')[-1].split('.')[0])}
        res.append(ds.copy() )

    return res


    