import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import argparse

import ppls1.imp.chp as chp

ap = argparse.ArgumentParser()
ap.add_argument("-p", "--path", required=True, type=str, help="path to binary checkpoint")
ap.add_argument("-e", "--export", required=True, type=str, help="export")
ap.add_argument("-q", "--quantity", required=True, nargs='+', help="quantity")
ap.add_argument("-n", "--nbins", required=False, default=11, type=int, help="number of bins")
args = vars(ap.parse_args())

fullPathChp = args['path']
expName = args['export']
quantity = args['quantity']
if type(quantity) is list:
    numQuantities = len(quantity)
else:
    quantity = [quantity]
    numQuantities = 1
nbins = args['nbins']

mass = 1.0

print('Path to simulation:       '+fullPathChp)
print('Path/name of pdf:         '+expName)
print('Quantities to be plotted: '+str(quantity))
print('Number of bins:           '+str(nbins))    


dfChp = chp.imp_chp_bin_DF(fullPathChp)

#dfChp = dfChp[dfChp['cid'] == 2]

bins = np.linspace(np.floor(dfChp['ry'].min()),np.ceil(dfChp['ry'].max()),nbins)
binVolume = (bins[1]- bins[0])*(np.ceil(dfChp['rx'].max())-np.floor(dfChp['rx'].min()))*(np.ceil(dfChp['rz'].max())-np.floor(dfChp['rz'].min()))

groups = dfChp.groupby(pd.cut(dfChp.ry, bins))

dfPlot = pd.DataFrame(index=np.diff(bins)/2+bins[:-1])
dfPlot['numPrtl'] = groups.ry.count()
dfPlot['rho_all'] = dfPlot['numPrtl']/binVolume
dfPlot['Tx_all'] = groups.vx.var()*mass
dfPlot['Ty_all'] = groups.vy.var()*mass
dfPlot['Tz_all'] = groups.vz.var()*mass
dfPlot['T_all'] = (dfPlot['Tx_all']+dfPlot['Ty_all']+dfPlot['Tz_all'])/3.0

if groups.ry.count().sum() != len(dfChp): print('Warning: groups.ry.count().sum() != len(dfChp)')

for idxQuant in range(numQuantities):
    if quantity[idxQuant] not in list(dfPlot.columns):
        print(f'Specified quantity {quantity[idxQuant]} not supported!')
        quantity[idxQuant] = 'numPrtl'

if numQuantities == 1:
    fig = plt.figure()
    (dfPlot[quantity[0]]).plot()
    plt.grid()
    plt.ylabel(quantity[0])
    plt.xlabel('z Coordinate')
    plt.xlim([np.floor(dfChp['ry'].min()),np.ceil(dfChp['ry'].max())])
else:  
    fig, axes = plt.subplots(numQuantities,1)
    #fig.suptitle(fullPathSim.split('/')[-2])
    #fig.tight_layout(pad=4.0)
    for idxQuant in range(numQuantities):
        (dfPlot[quantity[idxQuant]]).plot(ax=axes[idxQuant])
        axes[idxQuant].set_ylabel(quantity[idxQuant])
        axes[idxQuant].grid()
        axes[idxQuant].set_xlim([np.floor(dfChp['ry'].min()),np.ceil(dfChp['ry'].max())])
        if idxQuant < numQuantities-1: axes[idxQuant].set_xlabel('')
    axes[-1].set_xlabel('z Coordinate')

fig.savefig(expName, format='pdf')

print('Done')
