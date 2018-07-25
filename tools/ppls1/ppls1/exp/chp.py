# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 15:44:15 2018

@author: mheinen
"""

from struct import *

def exp_chp_bin(fname, chp):
    numParticles=len(chp['pid'])
    print('numParticles=',numParticles)
    with open(fname, "wb") as f:
        ba=bytearray()
        for pi in range(numParticles):
            ba.extend(pack('<Q',chp['pid'][pi]))
            ba.extend(pack('<I',chp['cid'][pi]))
            ba.extend(pack('<d',chp['rx'][pi]))
            ba.extend(pack('<d',chp['ry'][pi]))
            ba.extend(pack('<d',chp['rz'][pi]))
            ba.extend(pack('<d',chp['vx'][pi]))
            ba.extend(pack('<d',chp['vy'][pi]))
            ba.extend(pack('<d',chp['vz'][pi]))
            ba.extend(pack('<d',chp['q0'][pi]))
            ba.extend(pack('<d',chp['q1'][pi]))
            ba.extend(pack('<d',chp['q2'][pi]))
            ba.extend(pack('<d',chp['q3'][pi]))
            ba.extend(pack('<d',chp['Dx'][pi]))
            ba.extend(pack('<d',chp['Dy'][pi]))
            ba.extend(pack('<d',chp['Dz'][pi]))
        f.write(ba)

def exp_mmpld(chp, sphereparams):