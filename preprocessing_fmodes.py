import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
import datetime
from constants import DATA_PATH_DIRECTORY, DATA_PATH_OUTPUT, DATASET_HEADER, MIN_M, MAX_M, MIN_N, MAX_N, MIN_L, MAX_L
from data_analyser import find_all_frequencies
from data_reading import filename_finder_ell, read_frequency_dataset, savedata
import multiprocessing as mp
import math

# CODE TO CREATE AND SAVE TIMESERIES OF V1F DATA FOR A GIVEN AZIMUTHAL DEGREE (m), AND RADIAL ORDER (n),
# AND HARMONIC DEGREE (ell). 
# Unprocessed v1f data can be found on the GONG Helioseimic data archives: http://gong2.nso.edu/
# TO RUN THIS CODE, SAVE A DIRECTORY OF V1F DATA IN THE DIREWCTORY SPECIFIED BY datapath.
# THE DIRECTORY SHOULD CONTAIN ONLY THE V1F DATA, NOTHING ELSE.
# THE CODE WILL THEN SAVE THE .npz file IN THE DIRECTORY SPECIFIED BY DATA_PATH_OUTPUT, WHICH WILL CONTAIN
# THE TIMESERIES OF NU, DNU, AND THE AVG FREQ, AND MODE ANGLE FOR ALL COMBINATIONS OF m, n, and ell.

# THIS CODE USES data_analyser.py, data_reading.py, constants.py
# V1F DATA FROM mrv1f080702d***.txt AND LATER IN TIME HAVE 120 LINES OF HEADER (vs 119 FOR PREVIOUS DATASETS)
# AND THEREFORE skiprows MUST BE MANUALLY CHANGED AS NEEDED TO SUIT THE DATA YOU HAVE DOWNLOADED.  

# Authored: Tishtrya Mehta, 2023, tishtrya.mehta@gmail.com


#### FUNCTIONS ####

#----------------------------------------------------
# MODE_ANGLE
# Calculates the mode angle for a given m and ell using the formula cited in
# Simoniello 2016,  doi:10.3847/0004-637X/828/1/41
#    Parameters
#    ----------
#    m : int
#        The azimuthal degree
#    ell : int
#        The harmonic degree

#    Returns
#    -------
#    angle : float, or np.nan, np.nan
#        The mode angle for the given m and ell. If the mode angle is not possible, returns np.nan
#
def mode_angle(m, ell):
    if np.abs(m)/math.sqrt(ell*(ell+1)) <= 1:
        return math.degrees(math.acos(np.abs(m)/math.sqrt(ell*(ell+1))))
    else:
        return np.nan
#----------------------------------------------------

#----------------------------------------------------

# ANALYSIS_ELL
# Calculates the avg nu, nu time series array, and dnu time series array, and mode angle for all combinations of
# m, ell, and n. Takes in only a given value of ell, and returns all possible results for values between 
# MIN_M and MAX_M, and MIN_N and MAX_N, as defined in constants.py

#    Parameters
#    ----------
#    ell : int
#        The harmonic degree (between 20 and 150)
#    m : int
#        The azimuthal degree (between -ell and ell, defined in constants.py)
#    n : int
#        The radial order (between 1 and 35, defined in constants.py)

#    Returns
#    -------
#    temp_m_array : array
#        The m values for the given ell
#    temp_n_array : array
#        The n values for the given ell
#    temp_ell_array : array
#        The ell values for the given ell
#    temp_avg_nu_array : array
#        The avg nu (frequency) value for for the combination of ell, m, n
#    temp_nu_array : array
#        The nu (frequency) time series for the combination of ell, m, n
#    temp_dnu_array : array
#        The dnu (frequency) time series for the combination of ell, m, n
#    temp_mode_angle_array : array
#        The mode angle for the combination of ell, m, n
#
#
def analysis_ell(ell):
    dataset_names = filename_finder_ell(datapath, ell)
    folder_path= "GONG_v1f_"+f"{ell:03d}"+"_23"
    # THIS IS THE DIRECTORY WHERE ALL THE V1F DATA SHOULD BE STORED
    # NOTHING ELSE SHOULD BE IN THIS DIRECTORY


    # THIS IS CORRECT FOR CYCLE 23
    # [:-5] for cycle 23 as the final 5 datasets are not in the correct format and have skiprows = 119
    # You may need to change this depending on which mrv1y***.txt files you have downloaded from the GONG archives
    # and find the index of the first dataset with skiprows = 120

    all_frequency_datasets = {}
    for dataset_name in dataset_names[:-5]:
        all_frequency_datasets[dataset_name] = (
            read_frequency_dataset(DATA_PATH_DIRECTORY, folder_path, dataset_name, DATASET_HEADER, 119))
    for dataset_name in dataset_names[-5:]:
        all_frequency_datasets[dataset_name] = (
            read_frequency_dataset(DATA_PATH_DIRECTORY, folder_path, dataset_name, DATASET_HEADER, 120))
        

    temp_m_array=[]
    temp_n_array=[]
    temp_ell_array=[]
    temp_nu_array=[]
    temp_dnu_array=[]
    temp_avg_nu_array=[]
    temp_mode_angle_array=[]

    for m in np.arange(MIN_M, MAX_M+1, 1):
        for n in np.arange(MIN_N, MAX_N+1):
            array_nu, array_dnu = find_all_frequencies(m, n, all_frequency_datasets)
            temp_avg_nu = np.nanmean(array_nu)
            temp_mode_angle = mode_angle(m, ell)
            temp_avg_nu_array.append(temp_avg_nu)
            temp_m_array.append(m)
            temp_n_array.append(n)
            temp_ell_array.append(ell)
            temp_nu_array.append(array_nu)
            temp_dnu_array.append(array_dnu)
            temp_mode_angle_array.append(temp_mode_angle)

            #print("m: ", m, "n: ", n, "ell: ", ell, ", avg nu: ", temp_avg_nu, " COMPLETED")
    
    
    return temp_m_array, temp_n_array, temp_ell_array, temp_avg_nu_array, temp_nu_array, temp_dnu_array, temp_mode_angle_array

#----------------------------------------------------



#---------------------INPUTS------------------------
datapath = DATA_PATH_DIRECTORY+"/GONG_v1f_raw_23/"
outputname='/v1f_output_cycle_23.npz'
#----------------------------------------------------

#---------------------MAIN CODE------------------------
start_time = datetime.datetime.now()

# Find the length of the datasets by using a placeholder ell value
# and then use this to create empty arrays to store the nu and dnu data 
placeholder_ell=1
placeholder_datasets = filename_finder_ell(datapath, placeholder_ell)
save_nu_array=np.empty([0, len(placeholder_datasets)]) 
save_dnu_array=np.empty([0, len(placeholder_datasets)])

# Initialise empty arrays
save_m_array=[]
save_n_array=[]
save_ell_array=[]
save_avg_nu_array=[]
save_angle_array=[]
ell_array= np.arange(MIN_L,MAX_L+1)

# Run the analysis for each ell value in parallel
pool=mp.Pool()                                                                        
full_output=pool.map(analysis_ell,ell_array)
pool.close()
pool.join()

# Save the output
save_m_array=np.array([i[0] for i in full_output])
save_n_array=np.array([i[1] for i in full_output])
save_ell_array=np.array([i[2] for i in full_output])
save_avg_nu_array=np.array([i[3] for i in full_output])
save_nu_array=np.array([i[4] for i in full_output])
save_dnu_array=np.array([i[5] for i in full_output])
save_angle_array=np.array([i[6] for i in full_output])
# Flatten the arrays
save_m_array=save_m_array.flatten()
save_n_array=save_n_array.flatten()
save_ell_array=save_ell_array.flatten()
save_avg_nu_array=save_avg_nu_array.flatten()
save_nu_array=save_nu_array.reshape(-1, save_nu_array.shape[-1])
save_dnu_array=save_dnu_array.reshape(-1, save_dnu_array.shape[-1])
save_angle_array=save_angle_array.flatten()

# Save the output
np.savez(DATA_PATH_OUTPUT+outputname, m_array = save_m_array, n_array = save_n_array, ell_array = save_ell_array, avg_nu= save_avg_nu_array, nu_array = save_nu_array, dnu_array = save_dnu_array, angle = save_angle_array)
end_time = datetime.datetime.now()
# Print time taken
print("Time taken: ", end_time - start_time)
