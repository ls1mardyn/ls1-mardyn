# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 11:29:05 2018

@author: mheinen
"""
import numpy as np
from struct import *

#def imp_vdf_bin(files, flist=False):
def imp_vdf_bin(fplist):
    res=[]
#    fplist=[files]
#    if flist==True:
#        fplist=np.genfromtxt(files,dtype=np.str,)
        
#    print(fplist)
#    print("len:",len(fplist))
        
    for fp in fplist:
        f=open(fp,'rb')
        try:
            numClasses=(unpack('<L',f.read(4)))[0]
            numBins=(unpack('<L',f.read(4)))[0]
            dt=np.dtype('<u8')
            data=np.fromfile(f,dtype=dt)
            data=data.reshape(numClasses,numBins,order='F')
            # file info
            fpath=fp.rpartition('/')[0]
            fname=fp.rpartition('/')[-1]
            fext=fp.rpartition('.')[-1]
            finfo={'path':fpath,'name':fname,'ext':fext}
            quant=fname.split('_')[-2]
            flux=fname.split('_')[-3]
            comp=fname.split('_')[-4]
            reg=fname.split('_')[-5]
            distr=fname.split('_')[-6]
            ds={'numBins':numBins,
                'numClasses':numClasses,
                'data':data,
                'finfo':finfo,
                'quant':quant,
                'flux':flux,
                'comp':comp,
                'reg':reg,
                'distr':distr}
            res.append(ds.copy() )
            
        finally:
            f.close()
    return res


def imp_vdf_txt(fplist):
    res=[]
#    fplist=[files]
#    if flist==True:
#        fplist=np.genfromtxt(files,dtype=np.str,)
        
#    print(fplist)
#    print("len:",len(fplist))
        
    for fp in fplist:
        classes=np.genfromtxt(fp,dtype='<f8',usecols=0,skip_header=1)
        numClasses=len(classes)
        tmp=np.genfromtxt(fp,dtype=np.str,skip_footer=numClasses)
        bins=np.genfromtxt(fp,dtype='<f8',usecols=range(1,len(tmp)),skip_footer=numClasses)
        numBins=len(bins)

        dt=np.dtype('<u8')
        data=np.genfromtxt(fp,dtype=dt,skip_header=1,usecols=range(1,numBins+1))
        # file info
        fpath=fp.rpartition('/')[0]
        fname=fp.rpartition('/')[-1]
        fext=fp.rpartition('.')[-1]
        finfo={'path':fpath,'name':fname,'ext':fext}
        quant=fname.split('_')[-2]
        flux=fname.split('_')[-3]
        comp=fname.split('_')[-4]
        reg=fname.split('_')[-5]
        distr=fname.split('_')[-6]
        ds={'numBins':numBins,
            'numClasses':numClasses,
            'bins':bins,
            'classes':classes,
            'data':data,
            'finfo':finfo,
            'quant':quant,
            'flux':flux,
            'comp':comp,
            'reg':reg,
            'distr':distr}
        res.append(ds.copy() )
            
    return res


    