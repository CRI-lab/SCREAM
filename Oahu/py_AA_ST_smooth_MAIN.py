# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:31:56 2023

@author: rum

Changed the way we truncate the matrices to be more pythonic. 
TO DO: fix regression script.
Input: Text files / folders with the named text files for the region. 
Output: Figures and folders with linear regressions

ADAPTED from Tiffany Andersos AA_ST_smooth_MAIN.m' matlab script. 
This script calculates rates and intercepts using ST
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
#from AA_ST_regression_Python import AAA_ST_regress
from AA_ST_regression_Python import *
from make_smooth_mat_python import *

# BEGINNING OF USER INPUT AREA

flag_df_one_add_one = True
# flag_df_one_add_one = False  # Use this for Kauai...

region = 'EOahu'

savefig = True
# savefig = False

# END OF USER INPUT AREA
# ----------------------------------------------------------------------
m2ft = 3.28084  # conversion factor for meters to feet

# make directory called '<region>_results' if it does not exist
if not os.path.isdir(f"{region}_results"):
    os.mkdir(f"{region}_results")

# make directory called 'Figs' in existing dir '<region>_results'
# if it does not exist
if not os.path.isdir(f"{region}_results\Figs"):
    os.mkdir(f"{region}_results\Figs")

# if it does not exist
if not os.path.isdir(f"{region}_results\Rates"):
    os.mkdir(f"{region}_results\Rates")

# Define directory name and full filename where section file is located
regiondirname = f"{os.getcwd()}\{region}_dat\\"
sectfilename = f"{region}_sections.txt"
sectfile = f"{regiondirname}{sectfilename}"
sects = np.loadtxt(sectfile, dtype=str)

results_all = {}

for sectname in sects:#[5:8]:
    print(f"Processing {sectname}...")
    # Define directory name and full filename where data file is located
    dirname = f"{regiondirname}{sectname}_dat\\"
    modelname = f"{sectname}Toedist.txt"
    datfile = f"{dirname}{modelname}"

    # Define full filename of truncation file (identifies hardened shorelines)
    truncname = f"{sectname}Truncation.txt"
    truncfile = f"{dirname}{truncname}"

    # Define full filename of boundary file (identifies alongshore breaks in continuity)
    boundname = f"{sectname}Boundary.txt"
    boundfile = f"{dirname}{boundname}"

    # Define directory name and full filename where veg line data file is located
    vegname = f"{sectname}Vegdist.txt"
    vegfile = f"{dirname}{vegname}"

    # ---------------------------
    # load data, veg, boundary, and truncation files
    data_raw = pd.read_csv(datfile, sep="\t", header=None).to_numpy()
    full = np.copy(data_raw)
    full_flat = np.copy(data_raw)
    
    trunc = pd.read_csv(truncfile, sep="\t", header=None).to_numpy().astype(np.int16)
    bounds = np.loadtxt(boundfile)
    veg_raw = np.loadtxt(vegfile)

    if data_raw.size == 0:
        raise ValueError("Data file is empty. Please check file format. Text MS-DOS works.")
    if bounds.size == 0:
        raise ValueError("Boundary file is empty. Please do not use 'NaN' - define entire boundary.")
    if np.isnan(bounds).any():
        raise ValueError("Please do not use 'NaN' - define one boundary instead.")
    if trunc.size == 0:
        raise ValueError("Truncation file is empty. If no truncation, please use 'NaN'")
    if veg_raw.size == 0:
        raise ValueError("Veg dist file is empty. Please check file format. Text MS-DOS works.")
    if (veg_raw[2:, 0] == data_raw[2:, 0]).min() == False: ## ABM edit. dumb way to fix it
        print((veg_raw[2:, 0] == data_raw[2:, 0]).min())
        raise ValueError("Different number of survey dates in vegdist and toedist file. Please check data.")
    if not (veg_raw[2:, 0] == data_raw[2:, 0]).all():
        raise ValueError("Different number of transects in vegdist and toedist file. Please check data.")
    if not np.allclose(veg_raw[0, 1:], data_raw[0, 1:]):
        raise ValueError("Dates in veg dist file do not match those in toedist. Please check data.")
    if not (veg_raw[2:, 0] <= data_raw[2:, 0]).all():
        raise ValueError("Transect numbers in veg dist file do not match those in toedist. Please check data.")

    # define time series segments, years
    t_data = data_raw[0, 1:]
          
    # put truncated values in matrix for plotting later. initialize matrices. 
    hardshore_data = np.full_like(full, np.nan)                  #All nans except col_end value, and Nan'd values (full)
    truncated_data = np.full_like(full, np.nan)                  #All nans except Nan'd values (full)
    hardshore_data_flat = np.full_like(full, np.nan)             #All nans except col_end value, and Nan'd values  (full_flat)
    hardshore_data_flat_truncated = np.full_like(full, np.nan)   #All nans except Nan'd values (full_flat)
        
    
    ## Check that truncation inds are not just NaN, if so skip that file. No truncation
    if np.all(trunc==0):
        #raise ValueError("Nothing to truncate")
        print(f"Nothing to truncate {sectname}... trunc file empty")
        # continue # Skip to next elt in the loop
    else:
        ## Convert truncation inds to python (from matlab, columns index would start at 1)
        trunc[:, 2] -= 1
        
        # Based on truncation file, remove duplicate data due to seawall 
        # construction or other hard structure by setting those data to NaN 
        # NOTE: This does not change the size of the matrix (doesn't remove 
        # any rows or columns), it only replaces the duplicate data with NaNs or flattens it to the value in [row,col_end]
        
        for trans_start, trans_end, col_start, col_end in trunc:
            # Note that trans end is inclusive, but col end is not
            
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
            full[row_start:row_end, 1:col_start] = np.nan #if assuming third column is always 1, this is same as data_trunc

            full_flat[row_start:row_end, col_end:] = full[row_start:row_end, col_end-1].reshape((row_end-row_start, 1)) 
            
            #record second pair of data matrices
            hardshore_data_flat[row_start:row_end, col_end-1:] = full_flat[row_start:row_end, col_end-1:]
            hardshore_data_flat_truncated[row_start:row_end, col_end:] = full_flat[row_start:row_end, col_end:]
            
        
        #slice off the row with the dates and the column with the transect numbers.     
        hardshore_data = hardshore_data[2:,1:]
        truncated_data = truncated_data[2:,1:]
        hardshore_data_flat = hardshore_data_flat[2:,1:]
        hardshore_data_flat_truncated = hardshore_data_flat_truncated[2:,1:]
    
    ### END TRUNCATION SHIT, START REGRESSION
    
    # perform ST regression
    ST_sect, ST_sect_figs, ST_subset = AAA_ST_regress(full, bounds, flag_df_one_add_one)
    ST_sect['truncated_data_m'] = truncated_data
    ST_sect['truncated_data_ft'] = truncated_data * m2ft
    ST_sect['hardshore_data_m'] = hardshore_data
    ST_sect['hardshore_data_ft'] = hardshore_data * m2ft
    ST_sect['hardshore_data_flat_m'] = hardshore_data_flat
    ST_sect['hardshore_data_flat_ft'] = hardshore_data_flat * m2ft
    ST_sect['hardshore_data_flat_truncated_m'] = hardshore_data_flat_truncated
    ST_sect['hardshore_data_flat_truncated_ft'] = hardshore_data_flat_truncated * m2ft

    # save rates
    ST_subset.to_csv(f"{region}_results\Rates\{sectname}_rates.csv")

    if savefig:
        # save figures of data and parameters
        plt.figure(ST_sect_figs['data_METERS'])
        plt.savefig(f"{region}_results\Figs\{sectname}_{int(t_data[0])}_{int(t_data[-1])}_data_METERS.png")

        plt.figure(ST_sect_figs['params_METERS'])
        plt.savefig(f"{region}_results\Figs\{sectname}_{int(t_data[0])}_{int(t_data[-1])}_params_METERS.png")

        plt.figure(ST_sect_figs['data_FEET'])
        plt.savefig(f"{region}_results\Figs\{sectname}_{int(t_data[0])}_{int(t_data[-1])}_data_FEET.png")

        plt.figure(ST_sect_figs['params_FEET'])
        plt.savefig(f"{region}_results\Figs\{sectname}_{int(t_data[0])}_{int(t_data[-1])}_params_FEET.png")

    plt.close(ST_sect_figs['data_METERS'])
    plt.close(ST_sect_figs['params_METERS'])
    plt.close(ST_sect_figs['data_FEET'])
    plt.close(ST_sect_figs['params_FEET'])

    results_area = {
        'name': sectname,
        'ST': ST_sect,
        'dates': t_data.tolist(),
        'data_raw': data_raw.tolist(),
        'veg_raw': veg_raw.tolist()
    }

    results_all[sectname] = results_area

    # save workspace
    np.save(f"{region}_results\workspace_{sectname}", results_area)

np.save(f"{region}_results\workspace_ALL", results_all)
np.save(f"{region}_results\AA_results_all", results_all)

print(f"Pau: {region} rates.")

