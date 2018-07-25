# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 17:17:16 2018

@author: mheinen
"""

import pickle

def save_obj(fpath, obj):
    with open(fpath, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(fpath):
    with open(fpath, 'rb') as f:
        return pickle.load(f)
        