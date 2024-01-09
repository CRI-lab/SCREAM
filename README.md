# SCREAM
Shoreline Change Regression Erosion/Accretion Modeling (SCREAM). 

Goal of this repo is to take Tiffany Anderson's Matlab code to calculate Linear Regression via Historical imagery to Python. 
Note that this script does indexing differently from the original script to use Python indexing practices. 

SCRIPT: 
***AA_ST_smooth_MAIN_Python.py***
- Takes user input area to direct to the correct directory of historical modeling areas used by Anderson 2017. 
- Initializes directories for outputs if they don't exist.
- Uses "sectfile" to iterate through the different places
- Loads in the raw data (distances), truncation flats, and vegetation file.
- Creates four different new matrices with variations of the raw data depending on where the truncation is needed.
- Truncation file shows in this order: (transect row start, transect row end, column start to KEEP data, column start to KEEP data (inclusive)). If there is no truncation needed, it is "NaN".

These variables are for storage of what was changed in the original data matrix: 
* hardshore_data = Matrix is all nans except col_end value + Nan'd values.
* truncated_data = All nans except Nan'd values only (not including col_end)
* hardshore_data_flat = takes the value in col_end and changes the rest of the values to that value. (i.e: if data_raw(row, col_end) = 40.5, then data_raw(row, col_end:n) = 40.5). ((Where n is the total number of columns in the matrix))
* hardshore_data_flat_truncated = All nans, except with the same "flattened" values as in hardshore_data_flat.
  
and finally:

* **data_trunc** = What's used in the regression code. Full table copy with numbers and dates but with the NaNs that were truncated off.

Within AA_ST_smooth_MAIN_Python.py,  **AA_ST_regression_Python.py** is called! 


SCRIPT: 
***AA_ST_regression_Python.py***
* [] Finish debugging regression code
* !Note -- weights for time (line 60) seems too large...
* 

SCRIPT: 
***make_smooth_mat_python.py***
ChatGPT converted this script from Matlab to python. Used in regression code.
