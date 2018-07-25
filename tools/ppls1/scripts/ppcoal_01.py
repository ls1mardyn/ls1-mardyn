# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 10:14:29 2018

@author: mheinen
"""
import sys, os
from struct import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage

script_path=os.path.abspath(os.path.join(os.path.dirname(__file__)))
package_path=os.path.abspath(os.path.join(script_path, ".."))
print("ppls1 package path: ",package_path)
print("script path:        ",script_path)
sys.path.append(package_path)
import ppls1.imp.coal as imp
import ppls1.pp.coal as pp
import ppls1.plot.coal as pl

dts=2  # timestep in fs

#%% sim params
XZ=100
Y=100
shellWidth=6.0
binWidth=3.0

#%% import data
files1="flist1.txt"
files2="flist2.txt"
fplist1=np.genfromtxt(files1,dtype=np.str)
fplist2=np.genfromtxt(files2,dtype=np.str)
data1_sdd=imp.imp_coal_bin(fplist1[:], flist=True)
data2_sdd=imp.imp_coal_bin(fplist2[:], flist=True)
#data1=imp.imp_coal_txt(files1, flist=True)
#data2=imp.imp_coal_txt(files2, flist=True)

#%% process data
numBins=data1_sdd[0]['numBins']
numShells=data1_sdd[0]['numShells']
#ds1=data1[2]['data']
#ds2=data2[2]['data']
data3_sdd=data1_sdd.copy()
numDatasets=len(data1_sdd)
print("number of datasets: ",numDatasets)
rho_liq=np.empty(numDatasets,dtype=dict)
rho_vap=np.empty(numDatasets,dtype=dict)
for x in range(len(data1_sdd)):
    data3_sdd[x]['data']=pp.cat_secud(data1_sdd[x]['data'],data2_sdd[x]['data'])/6.022140857e-4
    
x,y=pp.grid_coords(XZ,Y,shellWidth,binWidth,numShells,numBins)
x=(x-Y*0.5)*0.1
y=(y-XZ*0.5)*0.1
#X,Y=np.meshgrid(x,y)

font = {'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 10,
}

#fig, ax = plt.subplots(nrows=1, ncols=1)   # contour plot
#fig = plt.figure(1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#for fi in range(0,10):
for fi in range(len(data3_sdd)):
#for d in data3:
    print("fi=",fi)

    data_sdd=data3_sdd[fi]['data']
    #data=ndimage.filters.gaussian_filter(data,[2,1])  # filter data
    
    fig=plt.pcolormesh(x, y, data_sdd, vmin=0., vmax=0.02/6.022140857e-4, cmap='jet', antialiased=False)
    ax=plt.gca()

    # contour plot
#    fig, ax = plt.subplots(nrows=1, ncols=1)   # contour plot
#    xm,ym=x[:-1],y[:-1]
#    cs_sdd=ax.contour(xm,ym,data_sdd,[15],linewidths=0.5,linestyles='solid',colors='black',label='sdd')
    
    plt.title(r'2D-Density field', fontdict=font)
        
    ax.set_xlabel(r'Coordinate $z$ / nm')
    ax.set_ylabel(r'Radius $r$ / nm')
    #ax.vlines(x[::10],0,960,color='white',lw=0.2)
    #ax.hlines(y[::50],0,1456.7832,color='white',lw=0.2)
    #ax.set_xlim(0,960)
    #ax.set_ylim(0,960)
    xlim=[-Y*0.5*0.1,Y*0.5*0.1]
    ylim=[-XZ*0.5*0.1,XZ*0.5*0.1]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    xg=[-50,50]
    yg=[-50,50]
    gc='grey'  # grid color
    pl.grid(ax,xg,yg,c=gc,lw=0.5)
    ax.set_aspect(1)
       
    # elapsed time label
    fs=8
    lc='white'
#    lc='black'
    ts=data3_sdd[fi]['timestep']
    etime=ts*1e-6*dts  # elapsed time in ns
    s=r'$t={:3.2f}$~ns'.format(etime)
    xt=-4
    yt=-4
    ax.text(xt, yt, s, fontsize=fs, color=lc)
        
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider=make_axes_locatable(ax)
    cax=divider.append_axes("right",size="5%",pad=0.1)
    cb=plt.colorbar(fig, cax=cax)
    cb.set_label(r'Density $\rho(z,r)$ / mol/l')
#   cb.set_label(r'Density deviation $(\rho_{dd}-\rho_{kdd})/(\rho_{dd})~[\%]$')
    
    ff='png'
    prefix='case01'
    fname="{2}_{0:09d}.{1}".format(ts,ff,prefix)
#   fname="coal_kdd{0:09d}.png".format(ts)
    plt.savefig(fname, format=ff, dpi=300, orientation='portrait', papertype='letter')
    #plt.show()  # needed to refresh!!!
