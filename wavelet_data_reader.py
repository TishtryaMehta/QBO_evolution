import datetime
import numpy as np
from scipy.signal import savgol_filter


#Extracts the Frequency shift form input data by subtracting the mean
def frequency_shift_caluculator(data):
    temp_mean=np.nanmean(data)
    subtracted=data-temp_mean
    return subtracted


def wavlet_data_reader(path, filename):
    #Read in the data
    data_full=np.genfromtxt(path+'/'+filename, skip_header=1, delimiter=',', dtype=float)
    data=[item[0] for item in data_full]

    #Ensure data is even  
    if (len(data)% 2) != 0:
        data=np.copy(data[:-1])
    return data


def data_preprocessing(data, dt, starttime, savgol_wl, savgol_pol):
    #Obtain frequency shifts
    data_shift=frequency_shift_caluculator(data)
    # Define the envelope of the data
    env=savgol_filter(data_shift,len(data_shift)-savgol_wl,savgol_pol)
    #substract the envelope from the data
    data_shift_norm=data_shift-env

    #Set up the time axis
    t=[starttime + datetime.timedelta(days =dt*i) for i in range(len(data_shift_norm))]
    return t, data_shift_norm, data_shift, env
