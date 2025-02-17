# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 15:44:15 2018
Modified 2021

@author: mheinen/homes
"""

from struct import pack

#%% Export binary ls1 checkpoint (representation: dict of lists)
def exp_chp_bin_DL(fname, chp, append=False):
    '''
    Export binary ls1 checkpoint (representation: dict of lists)

    :param string fname: Path and name of target
    :param chp: Checkpoint to be exported
    :param bool append: Specify if data should be appended to existing file; Default: False
    '''
    
    if append:
        writeMode = 'ab'
    else:
        writeMode = 'wb'
    
    numParticles=len(chp['pid'])
    
    with open(fname, writeMode) as f:
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

#%% Export binary ls1 checkpoint (representation: list of dict)
def exp_chp_bin_LD(fname, chp, append=False):
    '''
    Export binary ls1 checkpoint (representation: list of dict)

    :param string fname: Path and name of target
    :param chp: Checkpoint to be exported
    :param bool append: Specify if data should be appended to existing file; Default: False
    '''
    
    if append:
        writeMode = 'ab'
    else:
        writeMode = 'wb'
    
    with open(fname, writeMode) as f:
        ba=bytearray()
        for pi in chp:
#            li=list(pi.values())
#            print(*li)
            ba.extend(pack('<QIddddddddddddd',pi['pid'],pi['cid'],pi['rx'],pi['ry'],pi['rz'],pi['vx'],pi['vy'],pi['vz'],
                           pi['q0'],pi['q1'],pi['q2'],pi['q3'],pi['Dx'],pi['Dy'],pi['Dz'] ))
        f.write(ba)

#%% Export binary ls1 checkpoint (representation: data frame)
def exp_chp_bin_DF(fname, chp, append=False):
    '''
    Export binary ls1 checkpoint (representation: data frame)

    :param string fname: Path and name of target
    :param chp: Checkpoint to be exported
    :param bool append: Specify if data should be appended to existing file; Default: False
    '''
    
    if append:
        writeMode = 'ab'
    else:
        writeMode = 'wb'
    
    with open(fname, writeMode) as f:
        records = chp.to_records(index=False)
        records.tofile(f)
        f.close()

#%% Export ASCII ls1 checkpoint (representation: data frame)
def exp_chp_ascii_DF(fname, chp, meta, append=False):
    '''
    Export ASCII ls1 checkpoint (representation: data frame)

    :param string fname: Path and name of target
    :param chp: Checkpoint to be exported
    :param meta: Meta data incl. boxlengths to be exported
    :param bool append: Specify if data should be appended to existing file; Default: False
    '''
    
    if append:
        writeMode = 'ab'
    else:
        writeMode = 'wb'

    with open(fname, writeMode) as f:
        f.write("mardyn trunk 20120726\n")
        f.write("currentTime	0.0\n")
        f.write("Length	"+str(meta['boxlength'])+" "+str(meta['boxlength'])+" "+str(meta['boxlength'])+'\n')
        f.write("Temperature	"+str(5.0)+'\n')
        f.write("NumberOfComponents	1\n")
        f.write("1	0	0	0	0\n") ## Number of sites (LJ - Charge - Dipole - Quadrupole - Tersoff ??)
        f.write("0 0 0.0 1 1 0\n")   ## Site 1
        # f.write("0 0 0.4 1 1 0\n")   ## Site 2
        f.write("0.0 0.0 0.0\n")      ## Moment of inertia
        f.write("1e+10\n")
        f.write("NumberOfMolecules	"+str(len(chp))+'\n')
        f.write("MoleculeFormat	ICRVQD\n")
        
        for prtID in chp.index:
            f.write(str(prtID)+" 1 "+str(chp['x'][prtID])+" "+str(chp['y'][prtID])+" "+str(chp['z'][prtID])+"     "
                        +str(chp['vx'][prtID])+" "+str(chp['vy'][prtID])+" "+str(chp['vz'][prtID])+"     "
                        +str(chp['q1'][prtID])+" "+str(chp['q2'][prtID])+" "+str(chp['q3'][prtID])+" "+str(chp['q4'][prtID])+"     "
                        +str(chp['wx'][prtID])+" "+str(chp['wy'][prtID])+" "+str(chp['wz'][prtID])+'\n')
        f.close()

