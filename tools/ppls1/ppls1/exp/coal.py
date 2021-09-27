# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 11:53:04 2020

@author: mheinen
"""

import os
import json
from struct import *

# representation: dict of lists
def exp_coal_vtk2blob(outpath, data):
    for ds in data:
        # export json
        fname=ds['meta']['file_info']['name']; fprefix=os.path.splitext(fname)[0]
        fpath=f"{outpath}/{fprefix}.json"
        exp_json=ds['meta']
        exp_json['grid']['bins']['node_coords']=exp_json['grid']['bins']['node_coords'].tolist()
        exp_json['grid']['shells']['node_coords']=exp_json['grid']['shells']['node_coords'].tolist()
        with open(fpath,"w") as f:
            print(f"Exporting metadata to: {fpath}")
            json.dump(exp_json,f)

        fpath=f"{outpath}/{fprefix}.dat"
        with open(fpath, "wb") as f:
            print(f"Exporting data to: {fpath}")
            outdata=ds['data']
            ba=bytearray()
#            ba.extend('<Q',numBins)
#            ba.extend('<Q',numShells)
            ba.extend(outdata.tobytes(order='C'))
            f.write(ba)
