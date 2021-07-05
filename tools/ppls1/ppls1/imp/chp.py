# -*- coding: utf-8 -*-
'''
Created on Sat Jun  9 23:33:39 2018
Modified 2021

@author: mheinen/homes
'''

import os
import numpy as np
import pandas as pd
from struct import calcsize,unpack

#%% Import binary ls1 checkpoint (representation: dict of lists)
def imp_chp_bin_DL_legacy(fname, numParticles=99999999999, seek=0):
    '''
    LEGACY: Import binary ls1 checkpoint (representation: dict of lists)

    :param string fname: Path to binary checkpoint data
    :param int numParticles: Number of particles to be imported; Default: All
    :param int seek: Seek
    :return: chp: Dict of lists
    '''
    
    if numParticles == 99999999999:
        lengthInit = 0
    else:
        lengthInit = numParticles
    
    chp={}
    chp['pid']=np.zeros(lengthInit,dtype='<u8')
    chp['cid']=np.empty(lengthInit,dtype='<u4')
    chp['rx']=np.empty(lengthInit,dtype='<f8')
    chp['ry']=np.empty(lengthInit,dtype='<f8')
    chp['rz']=np.empty(lengthInit,dtype='<f8')
    chp['vx']=np.empty(lengthInit,dtype='<f8')
    chp['vy']=np.empty(lengthInit,dtype='<f8')
    chp['vz']=np.empty(lengthInit,dtype='<f8')
    chp['q0']=np.empty(lengthInit,dtype='<f8')
    chp['q1']=np.empty(lengthInit,dtype='<f8')
    chp['q2']=np.empty(lengthInit,dtype='<f8')
    chp['q3']=np.empty(lengthInit,dtype='<f8')
    chp['Dx']=np.empty(lengthInit,dtype='<f8')
    chp['Dy']=np.empty(lengthInit,dtype='<f8')
    chp['Dz']=np.empty(lengthInit,dtype='<f8')
    
    with open(fname, 'rb') as f:
        f.seek(seek)  # seek=0: go to beginning
        pi = 0
        while pi<numParticles:
            fmt='<QIddddddddddddd'
            buff_size=calcsize(fmt) # == 116 !!!
            buffer=f.read(buff_size)
            try:
                line=unpack(fmt,buffer)
                if numParticles == 99999999999:
                    chp['pid'] = np.append(chp['pid'],line[0])
                    chp['cid'] = np.append(chp['cid'],line[1])
                    chp['rx'] = np.append(chp['rx'],line[2])
                    chp['ry'] = np.append(chp['ry'],line[3])
                    chp['rz'] = np.append(chp['rz'],line[4])
                    chp['vx'] = np.append(chp['vx'],line[5])
                    chp['vy'] = np.append(chp['vy'],line[6])
                    chp['vz'] = np.append(chp['vz'],line[7])
                    chp['q0'] = np.append(chp['q0'],line[8])
                    chp['q1'] = np.append(chp['q1'],line[9])
                    chp['q2'] = np.append(chp['q2'],line[10])
                    chp['q3'] = np.append(chp['q3'],line[11])
                    chp['Dx'] = np.append(chp['Dx'],line[12])
                    chp['Dy'] = np.append(chp['Dy'],line[13])
                    chp['Dz'] = np.append(chp['Dz'],line[14])
                else:
                    chp['pid'][pi] = line[0]
                    chp['cid'][pi] = line[1]
                    chp['rx'][pi] = line[2]
                    chp['ry'][pi] = line[3]
                    chp['rz'][pi] = line[4]
                    chp['vx'][pi] = line[5]
                    chp['vy'][pi] = line[6]
                    chp['vz'][pi] = line[7]
                    chp['q0'][pi] = line[8]
                    chp['q1'][pi] = line[9]
                    chp['q2'][pi] = line[10]
                    chp['q3'][pi] = line[11]
                    chp['Dx'][pi] = line[12]
                    chp['Dy'][pi] = line[13]
                    chp['Dz'][pi] = line[14]
                pi += 1
            except Exception:
                break
        chp['pid'] = chp['pid'].astype('<u8', copy=False)
        chp['cid'] = chp['cid'].astype('<u4', copy=False)
        chp['rx'] = chp['rx'].astype('<f8', copy=False)
        chp['ry'] = chp['ry'].astype('<f8', copy=False)
        chp['rz'] = chp['rz'].astype('<f8', copy=False)
        chp['vx'] = chp['vx'].astype('<f8', copy=False)
        chp['vy'] = chp['vy'].astype('<f8', copy=False)
        chp['vz'] = chp['vz'].astype('<f8', copy=False)
        chp['q0'] = chp['q0'].astype('<f8', copy=False)
        chp['q1'] = chp['q1'].astype('<f8', copy=False)
        chp['q2'] = chp['q2'].astype('<f8', copy=False)
        chp['q3'] = chp['q3'].astype('<f8', copy=False)
        chp['Dx'] = chp['Dx'].astype('<f8', copy=False)
        chp['Dy'] = chp['Dy'].astype('<f8', copy=False)
        chp['Dz'] = chp['Dz'].astype('<f8', copy=False)
        
        if pi != numParticles:
            print(f'Warning: Less particles available than specified with numParticles. Importing {pi} particles ...')
            for key in chp: chp[key] = chp[key][:pi]
        
        #print(f'Imported {pi} particles')
    return chp

#%% Import binary ls1 checkpoint (representation: dict of lists)
def imp_chp_bin_DL(fname, numParticles=99999999999, seek=0):
    '''
    Import binary ls1 checkpoint (representation: dict of lists)

    :param string fname: Path to binary checkpoint data
    :param int numParticles: Number of particles to be imported; Default: All
    :param int seek: Seek
    :return: chp: Dict of lists
    '''
    
    dt = np.dtype([('pid', 'u8'), ('cid', 'u4'), ('rx', '<f8'), ('ry', '<f8'), ('rz', '<f8'), 
                   ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'), ('q0', '<f8'), ('q1', '<f8'),
                   ('q2', '<f8'), ('q3', '<f8'), ('Dx', '<f8'), ('Dy', '<f8'), ('Dz', '<f8')])
    data = np.fromfile(fname, dtype=dt)
    
    chp = {}
    lengthInit = min(data.size,numParticles-seek)
    chp['pid']=np.zeros(lengthInit,dtype='<u8')
    chp['cid']=np.empty(lengthInit,dtype='<u4')
    chp['rx']=np.empty(lengthInit,dtype='<f8')
    chp['ry']=np.empty(lengthInit,dtype='<f8')
    chp['rz']=np.empty(lengthInit,dtype='<f8')
    chp['vx']=np.empty(lengthInit,dtype='<f8')
    chp['vy']=np.empty(lengthInit,dtype='<f8')
    chp['vz']=np.empty(lengthInit,dtype='<f8')
    chp['q0']=np.empty(lengthInit,dtype='<f8')
    chp['q1']=np.empty(lengthInit,dtype='<f8')
    chp['q2']=np.empty(lengthInit,dtype='<f8')
    chp['q3']=np.empty(lengthInit,dtype='<f8')
    chp['Dx']=np.empty(lengthInit,dtype='<f8')
    chp['Dy']=np.empty(lengthInit,dtype='<f8')
    chp['Dz']=np.empty(lengthInit,dtype='<f8')
    
    pi = 0
    for dat in data:
        chp['pid'][pi] = dat[0]
        chp['cid'][pi] = dat[1]
        chp['rx'][pi] = dat[2]
        chp['ry'][pi] = dat[3]
        chp['rz'][pi] = dat[4]
        chp['vx'][pi] = dat[5]
        chp['vy'][pi] = dat[6]
        chp['vz'][pi] = dat[7]
        chp['q0'][pi] = dat[8]
        chp['q1'][pi] = dat[9]
        chp['q2'][pi] = dat[10]
        chp['q3'][pi] = dat[11]
        chp['Dx'][pi] = dat[12]
        chp['Dy'][pi] = dat[13]
        chp['Dz'][pi] = dat[14]
        pi += 1
    
    #print(f'Imported {pi} particles')
    return chp[seek:numParticles+seek]

#%% Import binary ls1 checkpoint (representation: list of dict)
def imp_chp_bin_LD(fname, numParticles=99999999999, seek=0):
    '''
    Import binary ls1 checkpoint (representation: list of dict)

    :param string fname: Path to binary checkpoint data
    :param int numParticles: Number of particles to be imported; Default: All
    :param int seek: Seek
    :return: chp: List of dict
    '''
    
    chp=[]
    
    with open(fname, 'rb') as f:
        f.seek(seek)  # seek=0: go to beginning
        pi = 0
        while pi<numParticles:
            fmt='<QIddddddddddddd'
            buff_size=calcsize(fmt) # == 116 !!!
            buffer=f.read(buff_size)
            try:
                line=unpack(fmt,buffer)
                mol={}
                mol['pid']=line[0]
                mol['cid']=line[1]
                mol['rx']=line[2]
                mol['ry']=line[3]
                mol['rz']=line[4]
                mol['vx']=line[5]
                mol['vy']=line[6]
                mol['vz']=line[7]
                mol['q0']=line[8]
                mol['q1']=line[9]
                mol['q2']=line[10]
                mol['q3']=line[11]
                mol['Dx']=line[12]
                mol['Dy']=line[13]
                mol['Dz']=line[14]
                chp.append(mol)
                pi += 1
            except Exception:
                break
        #print(f'Imported {pi} particles')
    return chp

#%% Import binary ls1 checkpoint (representation: data frame)
def imp_chp_bin_DF(fname, numParticles=99999999999, seek=0):
    '''
    Import binary ls1 checkpoint (representation: data frame)

    :param string fname: Path to binary checkpoint data
    :param int numParticles: Number of particles to be imported; Default: All
    :param int seek: Seek
    :return: chp: List of dict
    '''
    
    dt = np.dtype([('pid', 'u8'), ('cid', 'u4'), ('rx', '<f8'), ('ry', '<f8'), ('rz', '<f8'), 
                   ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'), ('q0', '<f8'), ('q1', '<f8'),
                   ('q2', '<f8'), ('q3', '<f8'), ('Dx', '<f8'), ('Dy', '<f8'), ('Dz', '<f8')])
    data = np.fromfile(fname, dtype=dt)
    return pd.DataFrame(data)[seek:numParticles+seek]

#%% Import ms2 restart file as data frame
def imp_ms2_rst_df(infile):
    '''
    Import ms2 restart file as data frame

    :param string infile: Path to ms2 restart file
    :return: df, meta: Data frame and meta data
    '''
    
    if infile == '':
        print('Filename empty!')
        return 0
    if os.path.exists(infile):
        with open(infile, 'r') as fileHandler:
            listOfLines = fileHandler.readlines()
    else:
        print('File does not exist!')
        return 0
    
    if listOfLines == 0: return 0
    meta = dict()
    meta['numPrtls'] = int(listOfLines[4])
    meta['boxlength'] = round((float(listOfLines[5]))**(1./3),8)
    meta['timestep'] = 0.0018236735
    df = pd.DataFrame(columns=['x','y','z','q1','q2','q3','q4','vx','vy','vz','wx','wy','wz'], index=range(meta['numPrtls']))
    # position
    listOfPos = listOfLines[13:13+meta['numPrtls']]
    for idx, val in enumerate(listOfPos):
        df['x'][idx] = meta['boxlength']*(0.5+float(val.split(';')[0]))
        df['y'][idx] = meta['boxlength']*(0.5+float(val.split(';')[1]))
        df['z'][idx] = meta['boxlength']*(0.5+float(val.split(';')[2]))
    # velos
    listOfVelo = listOfLines[13+meta['numPrtls']:13+2*meta['numPrtls']]
    for idx, val in enumerate(listOfVelo):
        df['vx'][idx] = (meta['boxlength']/meta['timestep'])*float(val.split(';')[0])
        df['vy'][idx] = (meta['boxlength']/meta['timestep'])*float(val.split(';')[1])
        df['vz'][idx] = (meta['boxlength']/meta['timestep'])*float(val.split(';')[2])
    # find position of quaternions
    meta['isElongated'] = False
    for idx, val in enumerate(listOfLines[13+2*meta['numPrtls']:]):
        if len(val.split(';')) == 4:
            meta['isElongated'] = True
            startQuad = 13+idx+2*meta['numPrtls']
            listOfQuad = listOfLines[startQuad:startQuad+meta['numPrtls']]
            for idx2, val2 in enumerate(listOfQuad):
                df['q1'][idx2] = float(val2.split(';')[0])
                df['q2'][idx2] = float(val2.split(';')[1])
                df['q3'][idx2] = float(val2.split(';')[2])
                df['q4'][idx2] = float(val2.split(';')[3])
            break
    # find position of angular velos
    if meta['isElongated']:
        for idx, val in enumerate(listOfLines[startQuad:]):
            if len(val.split(';')) == 3:
                meta['isElongated'] = True
                startW = startQuad+idx
                listOfW = listOfLines[startW:startW+meta['numPrtls']]
                for idx2, val2 in enumerate(listOfW):
                    df['wx'][idx2] = float(val2.split(';')[0])
                    df['wy'][idx2] = float(val2.split(';')[1])
                    df['wz'][idx2] = float(val2.split(';')[2])
                break
    else:
        df['q1'] = 1; df['q2'] = df['q3'] = df['q4'] = df['wx'] = df['wy'] = df['wz'] = 0
    return df, meta

