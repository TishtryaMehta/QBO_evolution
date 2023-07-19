import numpy as np
import pandas as pd
from constants import DATA_PATH_OUTPUT
import matplotlib.pyplot as plt
import datetime
from wavelet_data_reader import frequency_shift_caluculator
from v1f_analyser import index_fn, nan_excluder, frequency_sorted_mode
import os.path
import warnings



# MAIN
##########################

cycle_number=24

if cycle_number ==24:
    starttime=datetime.datetime(2008,12,29)
else:
    starttime=datetime.datetime(1996,9,3)


"""

#Load and read in data
v1f_data = np.load(DATA_PATH_OUTPUT+'/v1f_output_full_pool_angle_'+str(cycle_number)+'.npz', allow_pickle = True)
m_array = v1f_data['m_array']
n_array = v1f_data['n_array']
ell_array = v1f_data['ell_array']
nu_array = v1f_data['nu_array']
dnu_array = v1f_data['dnu_array']
avg_nu=v1f_data['avg_nu']
mode_angle=v1f_data['angle']

#index_n=index_fn(m=147, n=35, l=148)
#print(m_array[index_n])
#print(n_array[index_n])
#print(ell_array[index_n])
#print('**********')

"""

###############################
# Plot of latitude bands 
# for different frequency bands
###############################

LATITUDE_DATA_PATH= os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "CWT_analysis_fmodes", "Output", 'Sorted_v1f_csv_'+str(cycle_number)+'_98')

color_map = {
    3500.0: 'green',
    2700.0: 'orange',
    1900.0: 'red'
    
}
fig, ax = plt.subplots()

v1f_latitude = pd.read_csv(LATITUDE_DATA_PATH+ '/composite_results_cycle_'+str(cycle_number)+'_98.csv', header=0)
v1f_latitude.columns=v1f_latitude.columns.str.strip()


#plt.plot([], [], '^',  color='k', alpha=1.0, label='Avgd. over all latitudes')
#plt.plot(v1f_latitude['Avg_theta'],v1f_latitude['max_pow_period'], 'o',  color='blue', alpha=1.0)

# Iterate over each freq_str value and plot scatter points with error bars

for Avg_freq, group in v1f_latitude.groupby('Avg_freq'):
    # Filter out data with 'lat_avg' = 45
    #group = group[group['lat_avg'] != 37.51]
    x = group['Avg_theta']
    y = group['max_pow_period']
    y_err_lower = y - group['min_pow']  # Calculate the error bar heights for the lower ends
    y_err_upper = group['max_pow'] - y  # Calculate the error bar heights for the upper ends

    # Set error bar properties: black color with alpha=0.25 and capsize
    #ax.errorbar(x, y, yerr=[y_err_lower, y_err_upper], fmt='o', label=freq_str, color=color_map[freq_str], alpha=0.75, capsize=3)
    
    plt.fill_between(x, group['min_pow'], group['max_pow'], alpha=0.15, color=color_map[Avg_freq])
    plt.plot(x, y,  color=color_map[Avg_freq], alpha=1.0)

    #plot x,y, with label of Avg_freq, color of color_map[Avg_freq], alpha=1.0
    plt.plot(x, y, 'o', color=color_map[Avg_freq], alpha=1.0,  label=str(Avg_freq)+r' $\mu$Hz')


# Set plot title and labels
plt.title('Period corresponding to maximal CWT power against latitude')
plt.legend()
#yrange from 0 to 1500
plt.ylim(100,2000)
plt.xlabel('Latitude (degrees)')
plt.ylabel(r'Period (Days)')
plt.savefig(DATA_PATH_OUTPUT+'/QBO_latitude_'+str(cycle_number)+'_98.png', dpi=200)








##################
#BOTH CYCLES


"""
cycle_number=2324
seg_num='seg1'


color_map = {
    3500.0: 'green',
    2700.0: 'orange',
    1900.0: 'red'
    
}
fig, ax = plt.subplots()

LATITUDE_DATA_PATH= os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "CWT_analysis_fmodes", "Output", 'Sorted_v1f_csv_'+str(cycle_number)+'_98_'+seg_num)
v1f_latitude = pd.read_csv(LATITUDE_DATA_PATH+ '/composite_results_cycle_'+str(cycle_number)+'_98.csv', header=0)
v1f_latitude.columns=v1f_latitude.columns.str.strip()


#plt.plot([], [], '^',  color='k', alpha=1.0, label='Avgd. over all latitudes')
#plt.plot(v1f_latitude['Avg_theta'],v1f_latitude['max_pow_period'], 'o',  color='blue', alpha=1.0)

# Iterate over each freq_str value and plot scatter points with error bars

for Avg_freq, group in v1f_latitude.groupby('Avg_freq'):
    # Filter out data with 'lat_avg' = 45
    #group = group[group['lat_avg'] != 37.51]
    x = group['Avg_theta']
    y = group['max_pow_period']
    y_err_lower = y - group['min_pow']  # Calculate the error bar heights for the lower ends
    y_err_upper = group['max_pow'] - y  # Calculate the error bar heights for the upper ends

    # Set error bar properties: black color with alpha=0.25 and capsize
    #ax.errorbar(x, y, yerr=[y_err_lower, y_err_upper], fmt='o', label=freq_str, color=color_map[freq_str], alpha=0.75, capsize=3)
    
    plt.fill_between(x, group['min_pow'], group['max_pow'], alpha=0.15, color=color_map[Avg_freq])
    plt.plot(x, y,  color=color_map[Avg_freq], alpha=1.0)

    #plot x,y, with label of Avg_freq, color of color_map[Avg_freq], alpha=1.0
    plt.plot(x, y, 'o', color=color_map[Avg_freq], alpha=1.0,  label=str(Avg_freq)+r' $\mu$Hz')


seg_num='seg2'

LATITUDE_DATA_PATH= os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "CWT_analysis_fmodes", "Output", 'Sorted_v1f_csv_'+str(cycle_number)+'_98_'+seg_num)
v1f_latitude = pd.read_csv(LATITUDE_DATA_PATH+ '/composite_results_cycle_'+str(cycle_number)+'_98.csv', header=0)
v1f_latitude.columns=v1f_latitude.columns.str.strip()


for Avg_freq, group in v1f_latitude.groupby('Avg_freq'):
    # Filter out data with 'lat_avg' = 45
    #group = group[group['lat_avg'] != 37.51]
    x = group['Avg_theta']
    y = group['max_pow_period']
    y_err_lower = y - group['min_pow']  # Calculate the error bar heights for the lower ends
    y_err_upper = group['max_pow'] - y  # Calculate the error bar heights for the upper ends

    # Set error bar properties: black color with alpha=0.25 and capsize
    #ax.errorbar(x, y, yerr=[y_err_lower, y_err_upper], fmt='o', label=freq_str, color=color_map[freq_str], alpha=0.75, capsize=3)
    
    plt.fill_between(x, group['min_pow'], group['max_pow'], alpha=0.15, color=color_map[Avg_freq])
    plt.plot(x, y,  color=color_map[Avg_freq], alpha=1.0, linestyle='dashed')

    #plot x,y, with label of Avg_freq, color of color_map[Avg_freq], alpha=1.0
    plt.plot(x, y, '*', color=color_map[Avg_freq], alpha=1.0,  label=str(Avg_freq)+r' $\mu$Hz')




# Set plot title and labels
plt.title('Period corresponding to maximal CWT power against latitude')
plt.legend()
#yrange from 0 to 1500
plt.ylim(100,2000)
plt.xlabel('Latitude (degrees)')
plt.ylabel(r'Period (Days)')
plt.savefig(DATA_PATH_OUTPUT+'/QBO_latitude_'+str(cycle_number)+'_98.png', dpi=200)


"""
