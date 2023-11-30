# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:31:56 2023

@author: rum
"""

import pandas as pd
import numpy as np

import os


FILL_WITH_NANS = False

full_file_name = "EOahu_dat\Hauula_dat\HauulaToedist.txt"
trunc_file_name = "EOahu_dat\Hauula_dat\HauulaTruncation.txt"


## Load full file
full = pd.read_csv(full_file_name, sep="\t", header=None).to_numpy()

## Load truncation inds
trunc = pd.read_csv(trunc_file_name, sep="\t", header=None).to_numpy().astype(np.int16)

## Check that truncation inds are not NaN
if np.product(trunc.shape) == 1 and trunc.flatten()[0] is np.nan:
    
    raise ValueError("Nothing to truncate")

    #continue # Skip to next elt in the loop

## Convert truncation inds 

trunc[:, 2] -= 1

## Truncate relevant stuff

for trans_start, trans_end, col_start, col_end in trunc:
    # Note that trans end is inclusive, but col end is not
    
    print(trans_start, trans_end, col_start, col_end)
    
    row_start = np.where(full[:, 0] == trans_start)[0][0]
    row_end = np.where(full[:, 0] == trans_end)[0][0] + 1
    
    ## Since the first column (transect number) needs to be ignored
    col_start += 1 # We don't actually need column start at all!
    col_end += 1 
        
    ## Fill data
    
    if FILL_WITH_NANS:
        
        full[row_start:row_end, col_end:] = np.nan
        
    else:
        
        full[row_start:row_end, col_end:] = full[row_start:row_end, col_end-1].reshape((row_end-row_start, 1))
        
    print(full[row_start:row_end])
    
    
    
    

