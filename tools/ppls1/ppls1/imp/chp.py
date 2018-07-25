# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 23:33:39 2018

@author: mheinen
"""

import numpy as np
from struct import *

def imp_chp_bin(fname, numParticles):
    chp={}
    chp['pid']=np.empty(numParticles,dtype='<u8')
    chp['cid']=np.empty(numParticles,dtype='<u4')
    chp['rx']=np.empty(numParticles,dtype='<f8')
    chp['ry']=np.empty(numParticles,dtype='<f8')
    chp['rz']=np.empty(numParticles,dtype='<f8')
    chp['vx']=np.empty(numParticles,dtype='<f8')
    chp['vy']=np.empty(numParticles,dtype='<f8')
    chp['vz']=np.empty(numParticles,dtype='<f8')
    chp['q0']=np.empty(numParticles,dtype='<f8')
    chp['q1']=np.empty(numParticles,dtype='<f8')
    chp['q2']=np.empty(numParticles,dtype='<f8')
    chp['q3']=np.empty(numParticles,dtype='<f8')
    chp['Dx']=np.empty(numParticles,dtype='<f8')
    chp['Dy']=np.empty(numParticles,dtype='<f8')
    chp['Dz']=np.empty(numParticles,dtype='<f8')
    
    with open(fname, "rb") as f:
        f.seek(0)  # Go to beginning
        for pi in range(numParticles):
            fmt='<QIddddddddddddd'
            buff_size=calcsize(fmt) # == 116 !!!
            buffer=f.read(buff_size)
            line=unpack(fmt,buffer)
            chp['pid'][pi]=line[0]
            chp['cid'][pi]=line[1]
            chp['rx'][pi]=line[2]
            chp['ry'][pi]=line[3]
            chp['rz'][pi]=line[4]
            chp['vx'][pi]=line[5]
            chp['vy'][pi]=line[6]
            chp['vz'][pi]=line[7]
            chp['q0'][pi]=line[8]
            chp['q1'][pi]=line[9]
            chp['q2'][pi]=line[10]
            chp['q3'][pi]=line[11]
            chp['Dx'][pi]=line[12]
            chp['Dy'][pi]=line[13]
            chp['Dz'][pi]=line[14]           
    return chp