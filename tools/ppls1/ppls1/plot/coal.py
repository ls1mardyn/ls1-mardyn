# -*- coding: utf-8 -*-
"""
Created on Thu May 10 14:46:09 2018

@author: mheinen
"""

import matplotlib.pyplot as plt

def grid(ax,x,y,c,lw):
#    ax=plt.gca()
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    ax.vlines(x,ymin,ymax,color=c,lw=lw)
    ax.hlines(y,xmin,xmax,color=c,lw=lw)