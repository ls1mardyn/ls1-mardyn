# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 11:29:05 2018

@author: mheinen
"""
import numpy as np
import itertools
from operator import itemgetter

def vcat_pn(p,n):
    return np.vstack( (np.flipud(n),p) )

def hcat_pn(p,n):
    return np.hstack( (np.flipud(n),p) )

def classes(maxVal=None,numClasses=None):
    classWidth=maxVal/numClasses
    centers=np.linspace(0.5*classWidth,maxVal-0.5*classWidth,numClasses)
    centers_pn=np.hstack( (np.flipud(centers*-1),centers) ).tolist()
    res={'maxVal':maxVal,
         'numClasses':numClasses,
         'classWidth':classWidth,
         'centers':centers,
         'centers_pn':centers_pn}
    return res
    
def bins(Lz=None,numBins=None):
    binWidth=Lz/numBins
    centers=np.linspace(0.5*binWidth,Lz-0.5*binWidth,numBins)
    res={'Lz':Lz,
         'numBins':numBins,
         'binWidth':binWidth,
         'centers':centers}
    return res

def group_data(data):
    sdata=sorted(data, key=itemgetter('distr'))
    groups = {key:list(group) for key, group in itertools.groupby(sdata, key=lambda x:x['distr'])}
    
    for k in groups.keys():
        sdata = sorted(groups[k], key=itemgetter('reg'))
        groups[k] = {key:list(group) for key, group in itertools.groupby(sdata, key=lambda x:x['reg'])}
    
    for k1 in groups.keys():
        for k2 in groups[k1].keys():
            sdata = sorted(groups[k1][k2], key=itemgetter('comp'))
            groups[k1][k2] = {key:list(group) for key, group in itertools.groupby(sdata, key=lambda x:x['comp'])}
            
    for k1 in groups.keys():
        for k2 in groups[k1].keys():
            for k3 in groups[k1][k2].keys():
                sdata = sorted(groups[k1][k2][k3], key=itemgetter('flux'))
                groups[k1][k2][k3] = {key:list(group) for key, group in itertools.groupby(sdata, key=lambda x:x['flux'])}
    
    for k1 in groups.keys():
        for k2 in groups[k1].keys():
            for k3 in groups[k1][k2].keys():
                for k4 in groups[k1][k2][k3].keys():
                    sdata = sorted(groups[k1][k2][k3][k4], key=itemgetter('quant'))
                    groups[k1][k2][k3][k4] = {key:list(group) for key, group in itertools.groupby(sdata, key=lambda x:x['quant'])}    
    
    return groups
    