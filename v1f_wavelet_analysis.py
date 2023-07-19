import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from waveletFunctions import wavelet, wave_signif
import pylab
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import datetime
from matplotlib.dates import DateFormatter
from skimage.morphology import skeletonize, skeletonize_3d
from scipy.interpolate import interp1d
from scipy import ndimage
from csv import writer
import os.path
from datetime import date
from numpy.fft import fft, ifft
import math
import glob
import os
from v1f_analyser import index_fn, nan_excluder, frequency_sorted_mode, v1f_data_reader
from constants import DATA_PATH_OUTPUT, CSV_Header, DATA_PATH, SAVE_PATH
from wavelet_spectrum_maker import *
from v1f_data_sorting import wavelet_data_reader, data_preprocessing
from wavelet_errorbars import errorbar_finder
from wavelet_spectrum_plotter import spectrum_plotter
from wavelet_data_saver import wavelet_data_saver


##############################################################################################################
#PERFORM ANALYSIS ON PRESELECTED LATITUDE/FREQUENCY BANDS
##############################################################################################################

#load in csv file
cycle_number=2324
CYCLE_PATH= '/Sorted_v1f_csv_'+str(cycle_number)+'_98/'
CYCLE_PATH_SAVE= '/Sorted_v1f_csv_'+str(cycle_number)+'_98_seg2/'

if cycle_number ==24:
    starttime=datetime.datetime(2008,12,29)
    delete_number=1 #delete the first data point because the data has to be odd/even
else:
    starttime=datetime.datetime(1996,9,3)
    delete_number=0

if cycle_number ==2324:
    savgol_wl=100
    savgol_pol=7
else:
    savgol_wl=31
    savgol_pol=5

dt=36
seg1_index_min,seg1_index_max = 136, -1
index=0 #where the data is saved
DATA_PATH =  os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "CWT_analysis_fmodes", "Output")
#---------------------------


datalist_composite=[]

def loop_analysis(lower, higher, lower_theta, upper_theta):
    title='Wavelet_'+str(lower)+'_'+str(higher)+'_'+str(lower_theta)+'_'+str(upper_theta)
    data = pd.read_csv(DATA_PATH + CYCLE_PATH + str(lower)+'-'+str(higher)+'_'+str(lower_theta)+'-'+str(upper_theta)+'.csv', header=None)
    data=np.copy(data.to_numpy())
    data=data[~np.isnan(data)]
    data=np.copy(data[delete_number:])
    t, data_shift_norm, data_shift, env= data_preprocessing(data, dt, starttime, savgol_wl, savgol_pol)
    #WAVELET ANALYSIS
    time_array, period, scale, coi, power, global_ws, sig95, glsig95, normsig, global_signif= wavelet_analysis(data_shift_norm, dt)
    #CREATE ERRORBARS
    errorbar_packed= errorbar_finder(t, data, power, period, sig95, seg1_index_min, seg1_index_max)
    #PLOT
    spectrum_plotter(DATA_PATH + CYCLE_PATH_SAVE, title, t, data_shift, env, data_shift_norm, power, period, sig95, coi[1:], global_signif, global_ws, errorbar_packed)
    #SAVE
    datalist= wavelet_data_saver(DATA_PATH + CYCLE_PATH_SAVE, title,errorbar_packed,t)
    #datalist_composite.append(datalist)
    return datalist

lower=[1500,2300,3100]
higher=[2300,3100,3900]
lower_theta=[0,7,15,22,30,37,45,52,60,67]
upper_theta=[15,22,30,37,45,52,60,67,75,82]


for i in range(len(lower)):
    for j in range(len(lower_theta)):
        temp_datalist=loop_analysis(lower[i], higher[i], lower_theta[j], upper_theta[j])
        avg_freq= (lower[i]+higher[i])/2.
        avg_theta= (lower_theta[j]+upper_theta[j])/2.
        temp_extra=[lower[i], higher[i],avg_freq, lower_theta[j], upper_theta[j],avg_theta]
        temp_extra.extend(temp_datalist)
        datalist_composite.append(temp_extra)


file_exists = os.path.isfile(DATA_PATH + CYCLE_PATH_SAVE +'composite_results_cycle_'+str(cycle_number)+'.csv')
with open(DATA_PATH + CYCLE_PATH_SAVE +'composite_results_cycle_'+str(cycle_number)+'_98.csv', 'a') as f_object:
    if not file_exists:
        writer_object = writer(f_object)
        writer_object.writerow([CSV_Header])

    writer_object = writer(f_object)
    for i in range(len(datalist_composite)):
        writer_object.writerow(datalist_composite[i])
    f_object.close()
print('CSV output saved to:', DATA_PATH + CYCLE_PATH_SAVE + 'composite_results_cycle_'+str(cycle_number)+'.csv')

"""


##############################################################################################################
#PERFORM ANALYSIS ON ONE LATITUDE/FREQUENCY BAND
##############################################################################################################

cycle_number=23
filename= 'v1f_output_full_pool_angle_'+str(cycle_number)+'.npz'
save_path_ext= '/Output_lat_bands_'+str(cycle_number)+'/'
m_array, n_array, ell_array, nu_array, dnu_array, avg_nu, mode_angle= v1f_data_reader(DATA_PATH_OUTPUT, filename)
datalist_composite=[]

if cycle_number ==24:
    starttime=datetime.datetime(2008,12,29)
else:
    starttime=datetime.datetime(1996,9,3)



def loop_analysis(lower, higher, lower_theta, upper_theta, nu_array, mode_angle, ell_array, avg_nu, m_array, n_array):
    title='Wavelet_'+str(lower)+'_'+str(higher)+'_'+str(lower_theta)+'_'+str(upper_theta)
    data= frequency_sorted_mode(lower,higher, nu_array, mode_angle, lower_theta, upper_theta,  ell_array, avg_nu, m_array, n_array)
    data=data[~np.isnan(data)]
    data=np.copy(data[1:])
    t, data_shift_norm, data_shift, env= data_preprocessing(data, dt, starttime, savgol_wl, savgol_pol)
    #WAVELET ANALYSIS
    time_array, period, scale, coi, power, global_ws, sig95, glsig95, normsig, global_signif= wavelet_analysis(data_shift_norm, dt)
    #CREATE ERRORBARS
    errorbar_packed= errorbar_finder(t, data, power, period, sig95, seg1_index_min, seg1_index_max)
    #PLOT
    spectrum_plotter(SAVE_PATH+save_path_ext, title, t, data_shift, env, data_shift_norm, power, period, sig95, coi, global_signif, global_ws, errorbar_packed)
    #SAVE
    datalist= wavelet_data_saver(SAVE_PATH+save_path_ext, title,errorbar_packed,t)
    #datalist_composite.append(datalist)
    return datalist



lower=[1500,2300,3100]
higher=[2300,3100,3900]
lower_theta=[0,15,30,45,60]
upper_theta=[15,30,45,60,75]


lower=[2200]
higher=[2225]
lower_theta=[40]
upper_theta=[45]


for i in range(len(lower)):
    for j in range(len(lower_theta)):
        temp_datalist=loop_analysis(lower[i], higher[i], lower_theta[j], upper_theta[j], nu_array, mode_angle, ell_array, avg_nu, m_array, n_array)
        avg_freq= (lower[i]+higher[i])/2.
        avg_theta= (lower_theta[j]+upper_theta[j])/2.
        temp_extra=[lower[i], higher[i],avg_freq, lower_theta[j], upper_theta[j],avg_theta]
        temp_extra.extend(temp_datalist)
        datalist_composite.append(temp_extra)


file_exists = os.path.isfile(SAVE_PATH+save_path_ext+'composite_results_cycle_'+str(cycle_number)+'.csv')
with open(SAVE_PATH+save_path_ext+'composite_results_cycle_'+str(cycle_number)+'_single.csv', 'a') as f_object:
    if not file_exists:
        writer_object = writer(f_object)
        writer_object.writerow([CSV_Header])

    writer_object = writer(f_object)
    for i in range(len(datalist_composite)):
        writer_object.writerow(datalist_composite[i])
    f_object.close()
print('CSV output saved to:', SAVE_PATH+save_path_ext, 'composite_results_cycle_'+str(cycle_number)+'_single.csv')

"""

##############################################################################################################
#PERFORM ANALYSIS ON ONE LATITUDE/FREQUENCY BAND
##############################################################################################################

"""
#---------------------------
lower=1300    #1500-3900, 3100-3900, 1500-2300, 2300-3100
higher=1800
lower_theta=60       #7.5, 22.5, 37.5, 52.5, 67.5
upper_theta=90
savgol_wl=31
savgol_pol=5
starttime=datetime.datetime(1996,9,3)#datetime.datetime(2008,12,29)
dt=36
title='Wavelet_'+str(lower)+'_'+str(higher)+'_'+str(lower_theta)+'_'+str(upper_theta)
#contour_chooser=0
seg1_index_min,seg1_index_max = 0,-1
index=0 #where the data is saved
DATA_PATH =  os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "CWT_analysis_fmodes", "Output_lat_bands_24")
#---------------------------


filename= 'v1f_output_full_pool_angle_2324_20_85.npz'
m_array, n_array, ell_array, nu_array, dnu_array, avg_nu, mode_angle= v1f_data_reader(DATA_PATH_OUTPUT, filename)

data= frequency_sorted_mode(lower,higher, nu_array, mode_angle, lower_theta, upper_theta, ell_array, avg_nu)

#print(data)

#READ IN DATA
#data=wavelet_data_reader(DATA_PATH, filename, index)
#remove nans from data
#print('Number of nans removed', np.sum(np.isnan(data)))
data=data[~np.isnan(data)]
data=np.copy(data[1:])
#print(data)
#PREPROCESSING
t, data_shift_norm, data_shift, env= data_preprocessing(data, dt, starttime, savgol_wl, savgol_pol)
#WAVELET ANALYSIS
time_array, period, scale, coi, power, global_ws, sig95, glsig95, normsig, global_signif= wavelet_analysis(data_shift_norm, dt)
#CREATE ERRORBARS
errorbar_packed= errorbar_finder(t, data, power, period, sig95, seg1_index_min, seg1_index_max)
#PLOT
spectrum_plotter(SAVE_PATH, title, t, data_shift, env, data_shift_norm, power, period, sig95, coi, global_signif, global_ws, errorbar_packed)
#SAVE
wavelet_data_saver(SAVE_PATH, title,errorbar_packed,t)

#xx=latitude_bands_from_theta(30,75,50,1,15)
"""