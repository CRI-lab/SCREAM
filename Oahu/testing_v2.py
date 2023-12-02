# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:31:56 2023

@author: rum
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# BEGINNING OF USER INPUT AREA

flag_df_one_add_one = True
# flag_df_one_add_one = False  # Use this for Kauai...

region = 'EOahu'

savefig = True
# savefig = False

# END OF USER INPUT AREA
# ----------------------------------------------------------------------


full_file_name = "EOahu_dat\Hauula_dat\HauulaToedist.txt"
trunc_file_name = "EOahu_dat\Hauula_dat\HauulaTruncation.txt"


## Load full file with data (Toedist)
full = pd.read_csv(full_file_name, sep="\t", header=None).to_numpy() 
full_flat = np.copy(full)

## Load truncation inds 
trunc = pd.read_csv(trunc_file_name, sep="\t", header=None).to_numpy().astype(np.int16)

# put truncated values in matrix for plotting later. initialize matrices. 
hardshore_data = np.full_like(full, np.nan) #All nans except col_end value, and Nan'd values (full)
truncated_data = np.full_like(full, np.nan) #All nans except Nan'd values (full)
hardshore_data_flat = np.full_like(full, np.nan) #All nans except col_end value, and Nan'd values  (full_flat)
hardshore_data_flat_truncated = np.full_like(full, np.nan) #All nans except Nan'd values (full_flat)
    

## Check that truncation inds are not just NaN, if so skip that file. No truncation
if np.product(trunc.shape) == 1 and trunc.flatten()[0] is np.nan:
    
    raise ValueError("Nothing to truncate")

    #continue # Skip to next elt in the loop

## Convert truncation inds to python (from matlab, columns index would start at 1)

trunc[:, 2] -= 1

#Truncate relevant stuff. Based on trunctaion file, remove duplicate data due to seawall 
#     construction or other hard structure by setting those data to NaN 
#     NOTE: This does not change the size of the matrix (doesn't remove 
#           any rows or columns), it only replaces the duplicate data with NaNs or flattens it to the value in [row,col_end]

for trans_start, trans_end, col_start, col_end in trunc:
    # Note that trans end is inclusive, but col end is not
    
    #print(trans_start, trans_end, col_start, col_end)
    
    row_start = np.where(full[:, 0] == trans_start)[0][0]  # find the row index that matches the trans_start
    row_end = np.where(full[:, 0] == trans_end)[0][0] + 1  # find the row index that matches the trans_end + 1 (because need inclusivity)
    
    ## Since the first column (transect number) needs to be ignored
    col_start += 1 # We don't actually need column start at all!? It seems trunc file ALWAYS has third column starting at index 1 (to keep)
    col_end_og = col_end
    col_end += 1 
    
    #keep_cols = range(col_start,col_end,1)
    ## Fill data
    
    #record data matrices 
    hardshore_data[row_start:row_end, col_end-1:] = full[row_start:row_end, col_end-1:]
    truncated_data[row_start:row_end, col_end:] = full[row_start:row_end, col_end:]

    #change the data table copies
    full[row_start:row_end, col_end:] = np.nan #if assuming third column is always 1, this is same as data_trunc
    
    full_flat[row_start:row_end, col_end:] = full[row_start:row_end, col_end-1].reshape((row_end-row_start, 1)) 
    
    #record second pair of data matrices
    hardshore_data_flat[row_start:row_end, col_end-1:] = full_flat[row_start:row_end, col_end-1:]
    hardshore_data_flat_truncated[row_start:row_end, col_end:] = full_flat[row_start:row_end, col_end:]
    

#slice off the row with the dates and the column with the transect numbers.     
hardshore_data = hardshore_data[2:,1:]
truncated_data = truncated_data[2:,1:]
hardshore_data_flat = hardshore_data_flat[2:,1:]
hardshore_data_flat_truncated = hardshore_data_flat_truncated[2:,1:]

