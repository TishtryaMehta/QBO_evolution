from v1f_data_sorting import wavelet_data_reader, data_preprocessing
from wavelet_spectrum_maker import wavelet_analysis_finder
from constants_G_M import SAVE_PATH, DATA_PATH, DATA_PATH_DIRECTORY
from wavelet_errorbars import errorbar_finder
from wavelet_spectrum_plotter import spectrum_plotter
from wavelet_data_saver import wavelet_data_saver
import datetime
import os
import math
import numpy as np

# WRITTEN BY T MEHTA, 2023
# CONTACT TISHTRYA.MEHTA@GMAIL.COM IF YOU HAVE ANY PROBLEMS :)

# THIS CODE USES THE TORRENCE & COMPO WAVELET CODE, IN `waveletFunctions.py`
# AND OUTPUTS THE SIGNIFICANT PERIODS AND ERRORS IN PERIODICITY FOR THE
# INPUT DATA

##############################################################
# VARIABLES TO BE CHANGED BY THE USER
##############################################################

"""
#    User chosen variables
#    ----------
#    filename : string
#        Name of the file containing the data (csv or txt format)
#    index : int
#        The index of where the data is saved in the input file
#    savgol_wl : int
#        The window length for the savgol filter
#    savgol_pol : int
#        The polynomial order for the savgol filter
#    starttime : datetime.datetime
#        The start time of the data
#    dt : int
#        The time step of the data in days
#    title : string
#        The title of the output plot and csv file
#    seg1_index_min : int
#        The minimum index of the data to be examined
#    seg1_index_max : int
#        The maximum index of the data to be examined
#    DATA_PATH : string
#        The path to where the input data is saved
#    SAVE_PATH : string
#        Where the output plot and csv file are saved


This code returns
#    ----------
#    1. A plot of:
        i. The input data
        ii. The normalised data by divinding by a savgol envelope
        iii. The continuous wavelet spectrum of the normalised data with the 98% confidence level
            with the most significant period highlighted with errorbars
        iv. The global wavelet spectrum of the normalised data with the 95% confidence level
        
#    2. A csv file containing [columns are in order]
        a. The time coordinate of when the input data reaches its maximal value
        b. The time coordinate of when the input data reaches its minimal value
        c. The period corresponding to the maximal power observed in the continuous wavelet spectrum
        d. The period corresponding to the lower bound of the 98% confidence level for the maximal power period described in c.
        e. The period corresponding to the upper bound of the 98% confidence level for the maximal power period described in c.
        f. The time coordinate of when the input data reaches its maximal value at the 98% confidence level
        g. The period at the time described in a.
        h. The period corresponding to the lower bound of the 98% confidence level for the maximal power period described in g.
        i. The period corresponding to the upper bound of the 98% confidence level for the maximal power period described in g.
#    """

#---------------------------
#Default values
filename="GONG_low_24.dat"
savgol_wl=10
savgol_pol=5
starttime=datetime.datetime(2008,12,29)#datetime.datetime(1996,9,3)
dt=36
title='GONG_Low_Frequencies_Cycle_24'
seg1_index_min,seg1_index_max = 0,-1
index=0
DATA_PATH =  DATA_PATH
SAVE_PATH =  SAVE_PATH
#---------------------------

##############################################################
# MAIN
##############################################################

#READ IN DATA
#data=wavelet_data_reader(DATA_PATH_DIRECTORY, filename, index)
data_full=np.genfromtxt(DATA_PATH_DIRECTORY+filename)
data=[item[index] for item in data_full]
data=np.copy(data[:-1])

#remove nans from data
print('Number of nans removed', np.sum(np.isnan(data)))
data=data[~np.isnan(data)]
#PREPROCESSING
t, data_shift_norm, data_shift, env= data_preprocessing(data, dt, starttime, savgol_wl, savgol_pol)
#WAVELET ANALYSIS
time_array, period, scale, coi, power, global_ws, sig95, glsig95, normsig, global_signif= wavelet_analysis_finder(data_shift_norm, dt)
#CREATE ERRORBARS
errorbar_packed= errorbar_finder(t, data, power, period, sig95, seg1_index_min, seg1_index_max)
#PLOT
spectrum_plotter(SAVE_PATH, title, t, data_shift, env, data_shift_norm, power, period, sig95, coi, global_signif, global_ws, errorbar_packed)
#SAVE
wavelet_data_saver(SAVE_PATH, title,errorbar_packed,t)




         
