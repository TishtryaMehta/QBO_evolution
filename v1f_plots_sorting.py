import numpy as np
import pandas as pd
from constants import DATA_PATH_OUTPUT
import matplotlib.pyplot as plt
import datetime
from v1f_data_sorting import index_fn, frequency_sorted_mode, frequency_shift_caluculator, frequency_sorted_chooser
import os.path


##########################
# MAIN
##########################

cycle_number=23

if cycle_number ==24:
    starttime=datetime.datetime(2008,12,29)
else:
    starttime=datetime.datetime(1996,9,3)


#Load and read in data
v1f_data = np.load(DATA_PATH_OUTPUT+'/test_v1f.npz', allow_pickle = True)
m_array = v1f_data['m_array']
n_array = v1f_data['n_array']
ell_array = v1f_data['ell_array']
nu_array = v1f_data['nu_array']
dnu_array = v1f_data['dnu_array']
avg_nu=v1f_data['avg_nu']
mode_angle=v1f_data['angle']

###################################################
# Example case
###################################################
lower_freq, upper_freq=3100, 3500
lower_theta, upper_theta =0, 90
x= frequency_sorted_chooser(lower_freq,upper_freq, nu_array, mode_angle, lower_theta, upper_theta, ell_array, avg_nu, m_array, n_array)
t_GONG=[starttime + datetime.timedelta(days =36*i) for i in range(len(x))]
np.savetxt(DATA_PATH_OUTPUT+'/test_v1f_save.csv', x, delimiter=',')
plt.plot(t_GONG, x-np.mean(x),color='black', label=str(lower_theta)+'-'+str(upper_theta))
plt.legend()
plt.title(str(lower_freq)+'-'+str(upper_freq)+ ' uHz')
plt.xlabel('Years')
plt.ylabel(r'Frequency Shift ($\mu$Hz)')
plt.savefig(DATA_PATH_OUTPUT+'/test_v1f_plot.png')
#This data can now be manipulated/ plotted


###################################################
# Code to plot all frequency shifts for all 
# latitude bands
###################################################
"""
lower_freq=3100
upper_freq=3900

plt.figure()
lower_theta, upper_theta = 0,15
x= frequency_sorted_chooser(lower_freq,upper_freq, nu_array, mode_angle, lower_theta, upper_theta, ell_array, avg_nu, m_array, n_array)
t_GONG=[datetime.datetime(1996,9,3) + datetime.timedelta(days =36*i) for i in range(len(x))]
np.savetxt(DATA_PATH_OUTPUT+SORTED_PATH+str(lower_freq)+'-'+str(upper_freq)+'_'+str(lower_theta)+'-'+str(upper_theta)+'.csv', x, delimiter=',')
plt.plot(t_GONG, x-np.mean(x),color='black', label=str(lower_theta)+'-'+str(upper_theta))


lower_theta, upper_theta = 15,30
x= frequency_sorted_chooser(lower_freq,upper_freq, nu_array, mode_angle, lower_theta, upper_theta, ell_array, avg_nu, m_array, n_array)
np.savetxt(DATA_PATH_OUTPUT+SORTED_PATH+str(lower_freq)+'-'+str(upper_freq)+'_'+str(lower_theta)+'-'+str(upper_theta)+'.csv', x, delimiter=',')
plt.plot(t_GONG, x-np.mean(x),color='blue', label=str(lower_theta)+'-'+str(upper_theta))

lower_theta, upper_theta = 30,45
x= frequency_sorted_chooser(lower_freq,upper_freq, nu_array, mode_angle, lower_theta, upper_theta, ell_array, avg_nu, m_array, n_array)
np.savetxt(DATA_PATH_OUTPUT+SORTED_PATH+str(lower_freq)+'-'+str(upper_freq)+'_'+str(lower_theta)+'-'+str(upper_theta)+'.csv', x, delimiter=',')
plt.plot(t_GONG, x-np.mean(x),color='green', label=str(lower_theta)+'-'+str(upper_theta))

lower_theta, upper_theta = 45, 60
x= frequency_sorted_chooser(lower_freq,upper_freq, nu_array, mode_angle, lower_theta, upper_theta, ell_array, avg_nu, m_array, n_array)
np.savetxt(DATA_PATH_OUTPUT+SORTED_PATH+str(lower_freq)+'-'+str(upper_freq)+'_'+str(lower_theta)+'-'+str(upper_theta)+'.csv', x, delimiter=',')
plt.plot(t_GONG, x-np.mean(x),color='orange', label=str(lower_theta)+'-'+str(upper_theta))


lower_theta, upper_theta = 60,75
x= frequency_sorted_chooser(lower_freq,upper_freq, nu_array, mode_angle, lower_theta, upper_theta, ell_array, avg_nu, m_array, n_array)
np.savetxt(DATA_PATH_OUTPUT+SORTED_PATH+str(lower_freq)+'-'+str(upper_freq)+'_'+str(lower_theta)+'-'+str(upper_theta)+'.csv', x, delimiter=',')
plt.plot(t_GONG, x-np.mean(x),color='red', label=str(lower_theta)+'-'+str(upper_theta))

plt.legend()
plt.title(str(lower_freq)+'-'+str(upper_freq)+ ' uHz')
plt.xlabel('Years')
plt.ylabel(r'Frequency Shift ($\mu$Hz)')
plt.savefig(DATA_PATH_OUTPUT+str(lower_freq)+'-'+str(upper_freq)+'.png')
plt.close()
"""


##########################
# OTHER PLOTS
##########################

"""
#Histogram by nu_angle
plt.figure()
plt.hist(avg_nu[~np.isnan(avg_nu)], bins=25, color='blue', alpha=0.5,edgecolor='black', linewidth=1.2)
plt.xlabel(r'Average frequency ($\mu$Hz)')
plt.ylabel('Number of modes')
plt.savefig(DATA_PATH_OUTPUT+'/histogram_avg_nu_'+str(cycle_number)+'_all_ell.png')


#Histogram by mode angle
#excluded modes which have nan values in either 
#avg_nu or mode_angle
mode_angle_avg_nu=mode_angle[~np.isnan(avg_nu)]
mode_angle_avg_nu=mode_angle_avg_nu[~np.isnan(mode_angle_avg_nu)]
plt.figure()
plt.hist(mode_angle_avg_nu, bins=100, color='red', alpha=0.5,edgecolor='black', linewidth=1.2)
plt.xticks([0,15,30,45,60,75,90])
plt.xlabel(r'Mode Angle $\theta$ (degrees)')
plt.ylabel('Number of modes')
plt.savefig(DATA_PATH_OUTPUT+'/histogram_mode_angle_'+str(cycle_number)+'_all_ell.png')



#Nu-ell plot
plt.figure()
plt.scatter(ell_array[~np.isnan(avg_nu)], avg_nu[~np.isnan(avg_nu)], c=m_array[~np.isnan(avg_nu)], cmap='viridis', s=2.5)
#set colorbar label to be 'Azimuthal degree m'
cbar = plt.colorbar()
#colorbar label on the right hand side
#rotate the colorbar label
cbar.ax.set_ylabel('Azimuthal degree m', rotation=270)
cbar.ax.yaxis.set_label_position('right')
plt.xlabel(r'Harmonic Degree $\ell$')
plt.ylabel(r'Average frequency ($\mu$Hz)')
plt.savefig(DATA_PATH_OUTPUT+'/avg_nu_vs_ell_'+str(cycle_number)+'_all_ell.png')
"""
