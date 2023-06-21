# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 11:29:05 2018

@author: mheinen
"""
import os
import numpy as np
import json
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

def vtk_txt_extract_timestep(fpath, num_digits=7):
    fpath_noext=os.path.splitext(fpath)[0]
    number_str=fpath_noext[-num_digits::]
    timestep=int(number_str)
    return timestep

def imp_coal_vtk_txt(files, flist=False):
    res=[]
    fplist=[files]
    if flist==True:
        fplist=np.genfromtxt(files,dtype=np.str,)
        
#    print(fplist)
#    print("len:",len(fplist))

    ln_number=0
    points_params = list()
    cells_params = list()
    rho_params = list()
    with open(fplist[0]) as search:
        print('Scaning VTK format from file:',fplist[0])
        for line in search:
            ln_number+=1
            if line[0:6] =="POINTS":
                ln_points=ln_number
                points_params=line.split()
                print('reading POINTS data from line',ln_points)
            if line[0:5] =="CELLS":
                ln_cells=ln_number
                cells_params=line.split()
                print('reading CELLS data from line',ln_cells)
            if line[0] =="u":
                ln_rho=ln_number
                rho_params=line.split()
                print('reading u (rho) data from line',ln_rho)
                break
        
    for fp in fplist:
        print('Importing file:',fp)
        data_points=np.genfromtxt(fp,dtype=np.dtype('<f8'),skip_header=ln_points,max_rows=int(points_params[1]))
        data_cells=np.genfromtxt(fp,dtype=np.dtype('<u8'),skip_header=ln_cells,max_rows=int(cells_params[1]))
        data_rho=np.genfromtxt(fp,dtype=np.dtype('<f8'),skip_header=ln_rho,max_rows=int(rho_params[2]))
        indices=data_cells[:,1]
        data_rho=data_rho[indices]
        data_rho=data_rho[:,2]
        # determine num shells
        shell_coords=np.unique(data_points[:,1])
        shell_node_coords=shell_coords[::2]
        numShells=shell_node_coords.shape[0]-1
        maxShellBound=shell_node_coords[-1]
        shellHeight=shell_node_coords[1]-shell_node_coords[0]
        # determine num shells
        bin_coords=np.unique(data_points[:,0])
        bin_node_coords=bin_coords[::2]
        numBins=bin_node_coords.shape[0]-1
        maxBinBound=bin_node_coords[-1]
        binWidth=bin_node_coords[1]-bin_node_coords[0]
        binInfo={'node_coords':bin_node_coords,'count':numBins,'max_bound':maxBinBound,'width':binWidth}
        shellInfo={'node_coords':shell_node_coords,'count':numShells,'max_bound':maxShellBound,'height':shellHeight}
        gridInfo={'bins':binInfo,'shells':shellInfo}
        # reshape data to matrix
        data_rho_mat=np.reshape(data_rho, (-1,numBins))
        
        fpath=fp.rpartition('/')[0]
        fname=fp.rpartition('/')[-1]
        fext=fp.rpartition('.')[-1]
        finfo={'path':fpath,'name':fname,'ext':fext}
        timestep=vtk_txt_extract_timestep(fp,7)
        meta={'grid':gridInfo,
              'file_info':finfo,
              'timestep':timestep}
        ds={'data':data_rho_mat,
            'meta':meta}
        res.append(ds.copy() )
    return res

def imp_coal_vtk_blob(files, flist=False):
    res=[]
    fplist=[files]
    if flist==True:
        fplist=np.genfromtxt(files,dtype=np.str,)
        
#    print(fplist)
#    print("len:",len(fplist))
    for fp in fplist:
        fp_noext=os.path.splitext(fp)[0]
        with open(fp_noext+'.json',"r") as f:
            imp_json=json.load(f)
            imp_json['grid']['bins']['node_coords']=np.array(imp_json['grid']['bins']['node_coords'])
            imp_json['grid']['shells']['node_coords']=np.array(imp_json['grid']['shells']['node_coords'])
        with open(fp_noext+'.dat',"r") as f:
            dt=np.dtype('<f8')
            data=np.fromfile(f,dtype=dt)
            numShells,numBins=imp_json['grid']['shells']['count'],imp_json['grid']['bins']['count']
            data=data.reshape(numShells,numBins)
        ds={'meta':imp_json,
            'data':data}
        res.append(ds.copy() )
    return res