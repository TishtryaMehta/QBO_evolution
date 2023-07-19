import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime
from scipy.signal import savgol_filter


##########################
# FUNCTIONS
##########################


def frequency_sorted_chooser(lower,higher, nu_array, mode_angle, lower_theta, upper_theta,ell_array, avg_nu, m_array, n_array):
    """
    Sorts data by frequency and theta bands and outputs the average frequency time series dat. Only 
    outputs the avg of modes with 100% fill (i.e. no nans).

    Parameters
    ----------
    lower : int
        The lower frequency bound   
    higher : int
        The upper frequency bound
    nu_array : array
        The nu (frequency) time series for the combination of ell, m, n
    mode_angle : array
        The mode angle for the combination of ell, m, n
    lower_theta : int
        The lower theta bound
    upper_theta : int
        The upper theta bound
    ell_array : array   
        The ell values 
    avg_nu : array
        The avg nu (frequency) value for for the combination of ell, m, n
    m_array : array
        The m values
    n_array : array
        The n values

    Returns
    -------
    temp_nu_array : array
        The average frequency time series data, sorted by frequency and theta bands
    """
    index_frequency=np.where((avg_nu>=lower) & (avg_nu<higher) & (mode_angle<upper_theta) & (mode_angle>=lower_theta))
    temp_nu_array=[]
    ii_array=[]

    for ii in np.arange(len(nu_array[index_frequency][:,:])):
        if np.sum(np.isnan(nu_array[index_frequency][ii,:])) ==0:
            ii_array.append(ii)
            
    temp_nu_array=nu_array[index_frequency][ii_array]
    avg_angle=np.mean(mode_angle[index_frequency][ii_array])
    print("******************************************************************************")
    print('Between', lower, '--', higher, 'uHz and ', lower_theta, '--', upper_theta, 'deg, we averaged over ',len(temp_nu_array), ' modes, with avg angle:', avg_angle)
    # Uncomment to print out the modes that are being averaged over
    #print('l=',ell_array[index_frequency][ii_array], 'm=', m_array[index_frequency][ii_array], 'n=', n_array[index_frequency][ii_array], 'angle=',mode_angle[index_frequency][ii_array])
    return np.nanmean(temp_nu_array, axis=0)




def index_fn(m,n,l):
    """
    Finds the index of the given m,n,l in the v1f data array
    Parameters
    ----------
    m : int
        The azimuthal degree
    n : int
        The radial order    
    l : int
        The harmonic degree (between 20 and 150)

    Returns
    -------
    index : int
        The index of the given m,n,l in the v1f data array
        Or an error message if the index is out of bounds
    """
    temp_index=((l-20) * 35 * 301) + ((m + 150) * 35) + (n - 1)
    if (temp_index<0) or (temp_index>1380084) or (l<20) or (l>150) or (n<1) or (n>35) or (m<-l) or (m>l):
        print('INDEX ERROR!')
        print('Restrict l to be between 20 and 150')
        print('Restrict n to be between 1 and 35')
        print('Restrict m to be between -l and +l')
        print('*****************************************')
        return np.nan
    else:
        return ((l-20) * 35 * 301) + ((m + 150) * 35) + (n - 1)


#Function that reads in the data from the v1f files
def v1f_data_reader(datapath, v1f_filename):
    """
    Reads in the data from the v1f files
    Parameters
    ----------
    datapath : str
        The path to the directory containing the dataset
    v1f_filename : str  
        The path to the dataset relative to the directory

    Returns
    -------
    m_array : array
        The m values for the given ell
    n_array : array
        The n values for the given ell
    ell_array : array
        The ell values for the given ell
    nu_array : array
        The nu (frequency) time series for the combination of ell, m, n
    dnu_array : array
        The dnu (frequency) time series for the combination of ell, m, n
    avg_nu : array  
        The avg nu (frequency) value for for the combination of ell, m, n   
    mode_angle : array  
        The mode angle for the combination of ell, m, n 
    """
    v1f_data = np.load(datapath+'/'+v1f_filename, allow_pickle = True)
    m_array = v1f_data['m_array']
    n_array = v1f_data['n_array']
    ell_array = v1f_data['ell_array']
    nu_array = v1f_data['nu_array']
    dnu_array = v1f_data['dnu_array']
    avg_nu=v1f_data['avg_nu']
    mode_angle=v1f_data['angle']
    return m_array, n_array, ell_array, nu_array, dnu_array, avg_nu, mode_angle




def frequency_shift_caluculator(data):
    """
    Extracts the Frequency shift form input data by subtracting the mean
    Parameters
    ----------
    data : array
        The input data

    Returns
    -------
    subtracted : array
        The frequency shift data
    """
    temp_mean=np.nanmean(data)
    subtracted=data-temp_mean
    return subtracted


def wavelet_data_reader(path, filename, index):
    """
    Reads in the data from the wavelet files
    Parameters
    ----------
    path : str
        The path to the directory containing the dataset
    filename : str
        The path to the dataset relative to the directory
    index : int
        The index of the given m,n,l in the v1f data array

    Returns
    -------
    data : array
        The input data
    """
    data_full=np.genfromtxt(path+'/'+filename, skip_header=0, delimiter='/n', dtype=float)
    data=np.copy(data_full)
    #Ensure data is even  
    if (len(data)% 2) != 0:
        data=np.copy(data[:-1])
    return data


def data_preprocessing(data, dt, starttime, savgol_wl, savgol_pol):
    """
    Preprocesses the data by removing the envelope
    Parameters
    ----------
    data : array
        The input data
    dt : float
        The time step
    starttime : datetime
        The start time of the data
    savgol_wl : int
        The savgol window length
    savgol_pol : int
        The savgol polynomial order

    Returns
    -------
    t : array
        The time axis
    data_shift_norm : array
        The preprocessed data
    data_shift : array
        The frequency shift data
    env : array
        The envelope of the data    
    """

    data_shift=frequency_shift_caluculator(data)
    env=savgol_filter(data_shift,len(data_shift)-savgol_wl,savgol_pol)
    data_shift_norm=data_shift-env
    t=[starttime + datetime.timedelta(days =dt*i) for i in range(len(data_shift_norm))]
    return t, data_shift_norm, data_shift, env


def save_v1f(save_path, lower_freq, upper_freq, lower_theta, upper_theta, num_nans):
    """
    Saves the data to CSV for a given frequency and theta bands
    Parameters
    ----------
    save_path : str
        The path to the directory to save the data
    lower_freq : int
        The lower frequency bound
    upper_freq : int
        The upper frequency bound
    lower_theta : int
        The lower theta bound
    upper_theta : int
        The upper theta bound
    num_nans : int  
        The number of nans to exclude from the data
    """

    modes=frequency_sorted_mode(lower_freq[i],upper_freq[i], nu_array, mode_angle, lower_theta, upper_theta,num_nans, ell_array)
    np.savetxt(save_path+'/v1f_{}-{}_{}-{}.csv'.format(lower_freq[i],upper_freq[i],lower_theta, upper_theta), modes, delimiter=',')
    modes=[]


def theta_band_finder(theta_val,l_val):
    """
    #Function that finds the m value for a given theta value
    Parameters
    ----------
    theta_val : int
        The theta value
    l_val : int
        The harmonic degree (between 20 and 150)

    Returns
    -------
    m_val : int
        The m value for the given theta value
    """
    m_val=np.sqrt(l_val*(l_val+1))*math.cos(math.radians(theta_val))
    return np.abs(m_val)


##########################
# DEPRICATED FUNCTIONS!
##########################

def data_in_m_band(m_min,m_max, n_min,n_max):
    #search for directory bad_data if it exists, then read in the data from there    
    path = os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "Data_processed", "Azimuthal", "Full_data")
    nan_limit=5
    df=pd.DataFrame()

    directory = ["025","050","075","100","125"]
    for d_temp in directory:
        for n_temp in np.arange(n_min,n_max+1):
            for m_temp in np.arange(m_min,m_max+1):
                #print(m_temp)
                for file in glob.iglob(path+"_"+d_temp+'/v1f_m_'+str(m_temp)+'_n_'+str(n_temp)+'*.csv',recursive=True):
                    data_temp=np.genfromtxt(file, skip_header=1, delimiter=',', dtype=float)                    
                    data_temp_2=[item[0] for item in data_temp]
                    #if nan_count is greater than 30 (i.e. more than 30% of the data is nan), then skip this file
                    if np.count_nonzero(np.isnan(data_temp_2))>nan_limit:
                        continue
                    else:
                        df_temp=pd.DataFrame(data_temp_2)
                        df=pd.concat([df,df_temp],axis=1)
            for m_temp in np.arange(-1*m_min,-1*m_max-1, -1):
                #print(m_temp)
                for file in glob.iglob(path+'/v1f_m_'+str(m_temp)+'_n_'+str(n_temp)+'*.csv',recursive=True):
                    data_temp=np.genfromtxt(file, skip_header=1, delimiter=',', dtype=float)
                    data_temp_2=[item[0] for item in data_temp]
                    if np.count_nonzero(np.isnan(data_temp_2))>nan_limit:
                        continue
                    else:
                        df_temp=pd.DataFrame(data_temp_2)
                        df=pd.concat([df,df_temp],axis=1)
    return df


def latitude_bands_from_theta(theta_min,theta_max,l_val,nmin,nmax):
    #Obtain the m values for the given theta values
    m_value_1=round(theta_band_finder(theta_min,l_val))
    m_value_2=round(theta_band_finder(theta_max,l_val))

    m_array=sorted([m_value_1,m_value_2])
    m_min=m_array[0]
    m_max=m_array[1]
    print(m_min,m_max, 'm_min, m_max')

    #Read in data for the given range of m values
    temp_data_frame=data_in_m_band(m_min,m_max,nmin,nmax)
    print(temp_data_frame,'temp_data_frame')
    #print(temp_data_frame,'temp_data_frame')
    #Average across each row in dataframe, ignoring NAN values
    temp_data_frame['mean']=temp_data_frame.mean(axis=1,skipna=True)
    #print(temp_data_frame['mean'],'mean')

    #Save data
    savepath="/home/space/phrsnz/Desktop/Academic/Programs/PhD/QBO_evolution/Data_processed/Azimuthal/Latitude_bands/"
    temp_data_frame['mean'].to_csv(savepath+str(theta_min)+'_'+str(theta_max)+'.csv',index=False,header=False)

    fig, (ax1, ax2) = plt.subplots(2)
    #plot the mean

    ax1.plot(temp_data_frame['mean'])
    #ax1.set_title('Latitude band: '+str(theta_min)+'-'+str(theta_max)+' degrees, NAN limit=10 ')
    ax1.set_title('all l, all n, theta: 30-75')
    #plot the frequency shift, by subtracting the mean
    ax2.plot(frequency_shift_caluculator(temp_data_frame['mean']))
    #ax2.set_title('Freq shift: Latitude band: '+str(theta_min)+'-'+str(theta_max)+' degrees')
    #save to savepath
    plt.savefig(savepath+'L_ALL_'+str(theta_min)+'_'+str(theta_max)+'.png')
    print("Computation complete")



#Excludes data with more than num_nans nan values
dat_array=[]
def nan_excluder(inputdata, num_nans):
    for dat in inputdata:
        if np.sum(np.isnan(dat))<num_nans:
            dat_array.append(dat)
    print('number of datasets with >6 nans:', len(dat_array), 'out of', len(inputdata))
    return dat_array


def frequency_sorted_mode(lower,higher, nu_array, mode_angle, lower_theta, upper_theta, num_nans, ell_array, avg_nu, m_array, n_array):
#Sorts data by frequency and theta bands and outputs the average frequency time series dat
    index_frequency=np.where((ell_array>=20) & (avg_nu>=lower) & (avg_nu<higher) & (mode_angle<upper_theta) & (mode_angle>=lower_theta) & (np.isnan(avg_nu)==False))
    print(index_frequency, 'index_frequency')
    avg_angle=np.nanmean(mode_angle[index_frequency])
    nan_excluded_data = nan_excluder(nu_array[index_frequency][:,:], num_nans)
    print(np.shape(nan_excluded_data))
    print("******************************************************************************")
    print('Between', lower, '--', higher, 'uHz and ', lower_theta, '--', upper_theta, 'deg, we averaged over ',num_nans, ' modes, with avg angle:', avg_angle)
    nu_array_avg=np.nanmean(nan_excluded_data, axis=0)
    return nu_array_avg[2:]