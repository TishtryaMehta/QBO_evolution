import os

#----------------------------------------------------
# CONSTANTS

# The minimum and maximum values for m, n, and ell, used in preprocessing_fmodes.py
MIN_M = -150
MAX_M = 150
MIN_N = 1
MAX_N = 35
MIN_L= 20
MAX_L= 150
#----------------------------------------------------


#----------------------------------------------------
# PATHS

# Where the input data is stored, hosting the mrv1fd***.txt files
DATA_PATH_DIRECTORY = os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "Data_Raw")

# Where the output .npz file created by preprocessing_fmodes.py is stored
DATA_PATH_OUTPUT =  os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Wavelet_QBO_files")

DATA_PATH =  os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Wavelet_QBO_files")
SAVE_PATH = os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Wavelet_QBO_files")

#----------------------------------------------------


#----------------------------------------------------
# HEADERS

# The header of the mrv1fd***.txt files
DATASET_HEADER = [
    "n",
    "l",
    "m",
    "nuguess",
    "fwguess",
    "nu",
    "dnu",
    "fwhm",
    "dpfwhm",
    "dmfwhm",
    "psamp",
    "dppsamp",
    "dmpsamp",
    "bkg0",
    "dbkg0",
    "bkg1",
    "dbkg1",
    "merit",
    "niter",
    "bad",
    "ierr",
    "fhbad",
]

# The header written in v1f_wavelet_analysis.py
CSV_Header=[
    "Low_freq",
    "High_freq",
    "Avg_freq",
    "Low_theta",
    "High_theta",
    "Avg_theta",
    "Solarmax",
    "Solarmin",
    "max_pow_period",
    "min_pow",
    "max_pow",
    "t_max_pow_period",
    "period_solar_max",
    "min_pow_solar_max",
    "max_pow_solar_max",
    "index_max"
]

#----------------------------------------------------