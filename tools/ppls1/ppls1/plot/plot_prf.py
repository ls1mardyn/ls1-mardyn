from matplotlib import pyplot as plt
import argparse

import ppls1.pp.pp_prf2df as pp_prf2df

flgArgs = 1

if flgArgs == 1:
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--path", required=True, type=str, help="path to profile data")
    ap.add_argument("-e", "--export", required=True, type=str, help="export")
    ap.add_argument("-q", "--quantity", required=True, nargs='+', help="quantity")
    ap.add_argument("-a", "--avg", action='store_true', help="quantity")
    ap.add_argument("-t", "--timestep", type=int, help="quantity")
    args = vars(ap.parse_args())

    fullPathSim = args['path']
    expName = args['export']
    quantity = args['quantity']
    if type(quantity) is list:
        numQuantities = len(quantity)
    else:
        quantity = [quantity]
        numQuantities = 1
    flgAverage = args['avg']
    timestep = args['timestep']
else:
    fullPathSim = 'PATHTOSIMULATION'
    expName = 'test.pdf'
    quantity = 'rho_all'
    numQuantities = 1
    flgAverage = True
    timestep = None
    

print('Path to simulation:       '+fullPathSim)
print('Path/name of pdf:         '+expName)
print('Quantities to be plotted: '+str(quantity))
print('flgAverage:               '+str(flgAverage))
print('Timestep:                 '+str(timestep))

if (timestep is None and not flgAverage):
    print('Neither timestep nor averaging set. Defaulting to averaging...')
    flgAverage = True
    timestep = None

if (timestep is not None and flgAverage):
    print('Both timestep and averaging set. Defaulting to averaging...')
    flgAverage = True
    timestep = None
    


df, md = pp_prf2df.prf2df(fullPathSim, flgExport=False, quietMode=True)

dfOne = df[df['cid']==0]

if flgAverage:
    dfAvg = dfOne.groupby('pos').mean()
else:
    dfAvg = dfOne[dfOne['timestep'] == timestep].copy()

dfAvg.drop(columns=['cid','timestep'],inplace=True)


if numQuantities == 1:
    fig = plt.figure()
    (dfAvg[quantity[0]]).plot()
    plt.grid()
    plt.ylabel(quantity[0])
else:  
    fig, axes = plt.subplots(numQuantities,1)
    #fig.suptitle(fullPathSim.split('/')[-2])
    #fig.tight_layout(pad=4.0)
    for idxQuant in range(numQuantities):
        (dfAvg[quantity[idxQuant]]).plot(ax=axes[idxQuant])
        axes[idxQuant].set_ylabel(quantity[idxQuant])
        axes[idxQuant].grid()
        if idxQuant < numQuantities-1: axes[idxQuant].set_xlabel('')
        
#plt.show()
fig.savefig(expName, format='pdf')

print('Done')
