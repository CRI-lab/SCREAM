import os
import numpy as np
import matplotlib.pyplot as plt

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

    # put truncated values in matrix for plotting later
    truncated_data = np.full_like(data_raw, np.nan)
    # put hard shoreline values in matrix for plotting later
    hardshore_data = np.full_like(data_raw, np.nan)
    hardshore_data_flat = np.full_like(data_raw, np.nan)
    hardshore_data_flat_truncated = np.full_like(data_raw, np.nan)

    # Based on truncation file, remove duplicate data due to seawall
    # construction or other hard structure by setting those data to NaN
    # NOTE: This does not change the size of the matrix (doesn't remove
    #       any rows or columns), it only replaces the duplicate data with NaNs
    if np.isnan(trunc).all():
        # do nothing - no transect data to truncate
        pass
    else:
        for ind in range(trunc.shape[0]):
            row_start, _ = np.where(data_trunc[:, 0] == trunc[ind])  #CHANGE LATER
            row_start = row_start[-1] #CHANGE LATER
            row_end, _ = np.where(data_trunc[:, 0] == trunc[ind, 1]) #CHANGE LATER
            row_end = row_end[-1] #CHANGE LATER
            cols_keep = slice(trunc[ind, 2], trunc[ind, 3] + 1)
            cols_tmp = np.arange(data_trunc.shape[1] - 1)
            cols_tmp[cols_keep] = np.nan
            cols_remove = cols_tmp[~np.isnan(cols_tmp)] + 1

            if np.any(np.sum(~np.isnan(data_raw[row_start:row_end, trunc[ind, 3] + 1:]), axis=1) > 1):
                hard_dat = data_raw[row_start:row_end, trunc[ind, 3] + 1:]
                hard_dat_nan = np.isnan(hard_dat)
                hardshore_data[row_start:row_end, trunc[ind, 3] + 1:] = hard_dat
                hard_dat_flat = np.tile(hard_dat[:, 0], (1, hard_dat.shape[1]))
                hard_dat_flat[hard_dat_nan] = np.nan
                hardshore_data_flat[row_start:row_end, trunc[ind, 3] + 1:] = hard_dat_flat
                hardshore_data_flat_truncated[row_start:row_end, trunc[ind, 3] + 2:] = hard_dat_flat[:, 1:]

            truncated_data[row_start:row_end, cols_remove] = data_raw[row_start:row_end, cols_remove]
            data_trunc[row_start:row_end, cols_remove] = np.nan

    truncated_data = truncated_data[2:, 1:]
    hardshore_data = hardshore_data[2:, 1:]
    hardshore_data_flat = hardshore_data_flat[2:, 1:]
    hardshore_data_flat_truncated = hardshore_data_flat_truncated[2:, 1:]

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
