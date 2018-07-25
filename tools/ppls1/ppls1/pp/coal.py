# -*- coding: utf-8 -*-
"""
Created on Tue May  8 10:33:53 2018

@author: mheinen
"""
import numpy as np

def cat_secud(u,d):
    return np.vstack( (np.flipud(d),u) )

def grid_coords(XZ,Y,shellWidth,binWidth,numShells,numBins):
    shellWidth=XZ/np.floor(XZ/shellWidth)
    binWidth=Y/np.floor(Y/binWidth)
    ra=np.zeros(numShells)
    ra[0]=shellWidth
    ra[1]=np.sqrt(2)*ra[0]
    A1=ra[0]**2
    for s in range(2,numShells):
        ra[s]=np.sqrt(A1+ra[s-1]**2);
    
    dra=np.ones(numShells);
    dra[0]=ra[0];
    for s in range(1,numShells):
        dra[s]=ra[s]-ra[s-1];
    
    #x = np.linspace(-10, 10, numBins+1)
    #y = np.linspace(-10, 10, numShells+1)
    x=np.ones(numBins+1)*binWidth
    x[0]=0
    x=np.cumsum(x)
    y=np.vstack( (np.flipud(dra),dra) )
    y=np.insert(np.cumsum(y),[0],[0])
    return x,y