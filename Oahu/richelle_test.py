# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 13:04:07 2023

@author: rum
"""
import os
import numpy as np
import matplotlib.pyplot as plt

#testing script edited by Richelle
#11/22 
#sits in Oahu folder in historical_analysis_example. 
#references into the EOahu folders. 

# BEGINNING OF USER INPUT AREA

flag_df_one_add_one = True
# flag_df_one_add_one = False  # Use this for Kauai...

region = 'EOahu'

savefig = True
# savefig = False

# END OF USER INPUT AREA
# ----------------------------------------------------------------------


def truncate_data(truncation_file, data_table):
    # Load truncation file and data table
    # truncation_matrix = np.loadtxt(truncation_file, dtype=int)
    # data_matrix = np.loadtxt(data_table)
    data_matrix = data_table
    truncation_matrix = truncation_file;

    # Initialize the result matrix with None
    result_matrix = np.full_like(data_matrix, fill_value=None, dtype=data_matrix.dtype)

    # Initialize variables to keep track of changes
    modified_rows = set()
    modified_cols = set()

    # Iterate through each truncation row
    for start_row, end_row, start_col, end_col in truncation_matrix:
        # Ensure indices are within the data_matrix shape
        start_row = max(0, start_row - 1)
        end_row = min(data_matrix.shape[0], end_row)
        start_col = max(0, start_col - 1)
        end_col = min(data_matrix.shape[1], end_col)

        # Update the result_matrix with the selected data
        result_matrix[start_row:end_row, start_col:end_col] = data_matrix[start_row:end_row, start_col:end_col]

        # Update modified rows and columns
        modified_rows.update(range(start_row, end_row))
        modified_cols.update(range(start_col, end_col))

    return result_matrix, list(modified_rows), list(modified_cols)


m2ft = 3.28084  # conversion factor for meters to feet

# make directory called '<region>_results' if it does not exist
if not os.path.isdir(f"{region}_results"):
    os.mkdir(f"{region}_results")

# make directory called 'Figs' in existing dir '<region>_results'
# if it does not exist
if not os.path.isdir(f"{region}_results\Figs"):
    os.mkdir(f"{region}_results\Figs")

# Define directory name and full filename where section file is located
regiondirname = f"{os.getcwd()}\{region}_dat\\"
sectfilename = f"{region}_sections.txt"
sectfile = f"{regiondirname}{sectfilename}"
sects = np.loadtxt(sectfile, dtype=str)

results_all = {}

for sectname in sects:
    print(f"Processing {sectname}...")

    # Define directory name and full filename where data file is located
    dirname = f"{regiondirname}{sectname}_dat\\"
    modelname = f"{sectname}toedist.txt"
    datfile = f"{dirname}{modelname}"

    # Define full filename of truncation file (identifies hardened shorelines)
    truncname = f"{sectname}truncation.txt"
    truncfile = f"{dirname}{truncname}"

    # Define full filename of boundary file (identifies alongshore breaks in continuity)
    boundname = f"{sectname}boundary.txt"
    boundfile = f"{dirname}{boundname}"

    # Define directory name and full filename where veg line data file is located
    vegname = f"{sectname}vegdist.txt"
    vegfile = f"{dirname}{vegname}"

    # ---------------------------
    # load data, veg, boundary, and truncation files
    data_raw = np.loadtxt(datfile)
    trunc = np.loadtxt(truncfile)
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
    if not (veg_raw[0, 1:-1] == data_raw[0 , 1:-1]).all():
        raise ValueError("Different number of survey dates in vegdist and toedist file. Please check data.")
    if not (veg_raw[2:, 0] == data_raw[2:, 0]).all():
        raise ValueError("Different number of transects in vegdist and toedist file. Please check data.")
    if not np.allclose(veg_raw[0, 1:], data_raw[0, 1:]):
        raise ValueError("Dates in veg dist file do not match those in toedist. Please check data.")
    if not (veg_raw[2:, 0] <= data_raw[2:, 0]).all():
        raise ValueError("Transect numbers in veg dist file do not match those in toedist. Please check data.")

    # -----------------------------------------------
    # define time series segments
    t_data = data_raw[0, 1:]

    # Truncate full data series so that duplicate seawall positions are set to NaN
    data_trunc = np.copy(data_raw)

    # put truncated values in matrix for plotting later. initialize matrices. 
    truncated_data = np.full_like(data_raw, np.nan)
    # put hard shoreline values in matrix for plotting later
    hardshore_data = np.full_like(data_raw, np.nan)
    hardshore_data_flat = np.full_like(data_raw, np.nan)
    hardshore_data_flat_truncated = np.full_like(data_raw, np.nan)

    # Based on truncation file, remove duplicate data due to seawall
    # construction or other hard structure by setting those data to NaN
    # NOTE: This does not change the size of the matrix (doesn't remove
    #       any rows or columns), it only replaces the duplicate data with NaNs
    if np.isnan(trunc.flat[0]): # first element is NaN
        # do nothing - no transect data to truncate
        pass
    else:
        result, modified_rows, modified_cols = truncate_data(trunc, data_trunc)

# # Print result and modified indices
# print(result)
# print("Modified Rows:", modified_rows)
# print("Modified Columns:", modified_cols)

    # perform ST regression
    ST_sect, ST_sect_figs = AAA_ST_regress(data_trunc, bounds, flag_df_one_add_one)
    ST_sect['truncated_data_m'] = truncated_data
    ST_sect['truncated_data_ft'] = truncated_data * m2ft
    ST_sect['hardshore_data_m'] = hardshore_data
    ST_sect['hardshore_data_ft'] = hardshore_data * m2ft
    ST_sect['hardshore_data_flat_m'] = hardshore_data_flat
    ST_sect['hardshore_data_flat_ft'] = hardshore_data_flat * m2ft
    ST_sect['hardshore_data_flat_truncated_m'] = hardshore_data_flat_truncated
    ST_sect['hardshore_data_flat_truncated_ft'] = hardshore_data_flat_truncated * m2ft

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


