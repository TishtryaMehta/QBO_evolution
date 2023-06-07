from wavelet_data_reader import wavlet_data_reader, data_preprocessing
from wavelet_spectrum_maker import wavelet_analysis
from constants_wavelet import DATA_PATH, SAVE_PATH
from wavelet_errorbars import errorbar_finder
from spectrum_plotter import spectrum_plotter
from wavelet_data_saver import wavelet_data_saver
import datetime

# WRITTEN BY T MEHTA, 2023
# CONTACT TISHTRYA.MEHTA@GMAIL.COM IF YOU HAVE ANY PROBLEMS :)

# THIS CODE USES THE TORRENCE & COMPO WAVELET CODE, IN `waveletFunctions.py`
# AND OUTPUTS THE SIGNIFICANT PERIODS AND ERRORS IN PERIODICITY FOR THE
# QBO USING HELIOSESIMIC AND RADIO DATA
# CODE IS WRITTEN FOR PYTHON 3.6.9

##############################################################
# Read in data and sort according to datatype and duration
##############################################################

#---------------------------
m_val=-94
n_val=2
filename="v1f_m_"+str(m_val)+"_n_"+str(n_val)+"_interpolate.csv"
savgol_wl=31
savgol_pol=5
starttime=datetime.datetime(1996,9,3)
dt=36
title='m= '+str(m_val)+' n='+str(n_val)+'_cleanup'
contour_chooser=4
seg1_index_min,seg1_index_max = 0,-1
#---------------------------



#READ IN DATA
data=wavlet_data_reader(DATA_PATH, filename)
#PREPROCESSING
t, data_shift_norm, data_shift, env= data_preprocessing(data, dt, starttime, savgol_wl, savgol_pol)
#WAVELET ANALYSIS
time_array, period, scale, coi, power, global_ws, sig95, glsig95, normsig, global_signif= wavelet_analysis(data_shift_norm, dt)
#CREATE ERRORBARS
errorbar_packed= errorbar_finder(t, data, power, period, sig95, seg1_index_min, seg1_index_max,contour_chooser)
#PLOT
spectrum_plotter(SAVE_PATH, title, t, data_shift, env, data_shift_norm, power, period, sig95, coi, global_signif, global_ws, errorbar_packed)
#SAVE
wavelet_data_saver(SAVE_PATH, title,errorbar_packed,t)




         
