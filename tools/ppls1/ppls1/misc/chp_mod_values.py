import ppls1.imp.chp as imp
import ppls1.exp.chp as exp
import numpy as np

work_folder = 'PATHTOSIMULATION'
in_file_path = work_folder+'cp_binary_half.restart.dat'
out_file_path = work_folder+'cp_binary_full.restart.dat'

# Read in checkpoint
chp = imp.imp_chp_bin_LD(in_file_path)

partcomp = list()

# replicate particles
for par in chp:
    if par['cid'] == 1:
        par['cid'] = 2
    tmp = par.copy()
    tmp['ry'] = 569 - tmp['ry']
    par['ry'] = par['ry'] - 20 + 569
    partcomp.append(par)
    partcomp.append(tmp)
    
# refresh particle ids
for pi in range(len(partcomp)):
    partcomp[pi]['pid']=pi+1

a = list()
b = list()
for i in range(len(chp)):
    a.append(partcomp[i]['cid'])
    b.append(partcomp[i]['ry'])
print(set(a))

exp.exp_chp_bin_LD(out_file_path,partcomp)

#print(partcomp)
