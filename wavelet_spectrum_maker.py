from waveletFunctions import wavelet, wave_signif
import numpy as np



##############################################################
# DEFINE FUNCTIONS
##############################################################

#""""""""""""""""""""""""""""""""""""""""""""""""""""
# BROKEN POWER LAW
#""""""""""""""""""""""""""""""""""""""""""""""""""""

def broken_linear(x, a, b, c,d):
    return np.piecewise(x, [x < b, x >= b], [lambda x:-a*x, lambda x:(-(a*d)*b - d*x) + c])

#""""""""""""""""""""""""""""""""""""""""""""""""""""
# POWER SPECTRUM
# Calculates the power spectrum of the input data
#""""""""""""""""""""""""""""""""""""""""""""""""""""

def power_spectrum(x,dt):
    X = fft(x)
    n = np.arange(len(X))
    T = len(X)*dt
    freq = np.arange(len(X))/(len(X)*dt) 
    #alog freq and np.abs(X)
    alogfreq=np.log10(freq)
    alogX=np.log10(np.abs(X))
    x_data=1./alogfreq#(1./freq)
    y_data=alogX#(np.abs(X))
    #remove all infs and nans
    array1=np.where(np.isnan(x_data))
    array2=np.where(np.isnan(y_data))
    array3=np.where(np.isinf(x_data))
    array4=np.where(np.isinf(y_data))
    x_data2=np.delete(x_data,array3)
    y_data2=np.delete(y_data,array3)
    #remove first and final 3 pts correspinding to detrended region
    x_data2=np.copy(x_data2[3:])
    y_data2=np.copy(y_data2[3:])
    #Fitting
    popt, pcov = curve_fit(broken_linear, x_data2, y_data2, p0=[3,1,1,1])
    alpha=popt[0]*(-1)
    
    #Uncomment to plot
    #plt.figure()
    #plt.plot(1./alogfreq,alogX)
    #plt.xlabel('Period (Log$_{10}$)')
    #plt.ylabel(r'FFT Amplitude (Log$_{10}$)')
    #plt.plot(x_data2, (broken_linear(x_data2, *popt)), 'r-', label='fit')
    #plt.savefig("/home/space/phrsnz/Desktop/Academic/Programs/PhD/QBO_evolution/CWT_analysis_azimuthal/"+"TEST_DELETE.png", dpi=200)
    return freq, np.abs(X), alpha



#""""""""""""""""""""""""""""""""""""""""""""""""""""
# LAG1 FINDER
# Calculates the lag1 autocorrelation of the input data
#""""""""""""""""""""""""""""""""""""""""""""""""""""
def lag1_finder(data):
    x_dat = np.array(data)
    mean_dat = np.mean(x_dat)
    var_dat = np.var(x_dat)
    #Normalise
    ndata = x_dat - mean_dat
    #Autocorrelation
    acorr = np.correlate(ndata, ndata, 'full')[len(ndata)-1:] 
    acorr = acorr / var_dat / len(ndata)
    #return lag1 value
    return acorr[1]

##############################################################



#""""""""""""""""""""""""""""""""""""""""""""""""""""
# WAVELET ANALYSIS
#""""""""""""""""""""""""""""""""""""""""""""""""""""
def wavelet_analysis_finder(data, dt):
    #Calls by WAVELETFUNCTIONS.py (Obatined from GitHub by 'Evgeniya Predybaylo') 
    num=len(data)
    time_array=np.arange(num)*dt
    mother = 'MORLET' 
    # here we'll use a Morlet mother wavelet, 
    #this can be changed see notes in waveletFunctions.py
    variance = np.std(data, ddof=1) ** 2 # variance of the signal. 

    """
    Compute the wavelet transform
    -----------------------------
    Now that we have the signal we want to compute the wavelet transform.
    For this we'll use the `wavelet` function that takes arguments of the signal
    and the dt - time step/sampling time.
    It also takes keyword argments or optional inputs.
    Here, were basically only passing the ones we really care about, 
    and the rest will be set to the default values.
    You can pass things like a different mother wavelet, 
    a different mother wavelet parameter, 
    define the different scales, and pad to get up to correct number of scales etc. To
    see more about these, read the docstring above the `wavelet` function 
    in waveletFunctions. 
    Upon calling `wavelet` it returns the wavelet transform (wave), 
    the period, scale and the 
    cone-of-influence (coi). 
    """
    wave, period, scale, coi = wavelet(data, dt, mother="MORLET",  dj=1/30, pad=1)
    # compute transform
    power = (np.abs(wave)) ** 2
        # wavelet power is defined as absolute value of transform squared
    global_ws = (np.sum(power, axis=1) / num)
        # define global wavelet - which is time-average over all times

    """
    Compute the significance values
    -------------------------------
    To do these we'll use the `wave_signif` function. 
    Here you need to pass the signal or the variance of the signal,
    the sampling time, 
    and the scale (which is returned from above function).
    In this function you can also pass other keyword arguments including 
    the type of significance test you want to do (default chi-squared from 
    paper),a `lag1` parameter
    if you want to estimate red-noise, 
    a significance level (default 0.95/95%) and the 
    degrees of freedom.
    """
    #powerspecfreq,powerspecamp,alpha=power_spectrum(data,dt)
    #We determine the value of lag1 by running an autocorrelation in lag1_finder()
    lag1=lag1_finder(data)
    #print('Lag1 for this data is: ',lag1)
    signif = wave_signif(variance, dt=dt, sigtest=0, scale=scale, siglvl=0.98, lag1=lag1, mother=mother) 
    # calculate the significance 
    sig95 = signif[:, np.newaxis].dot(np.ones(num)[np.newaxis, :]) 
    # expand signif --> (J+1)x(N) array
    sig95 = power / sig95 
    # where ratio > 1, power is significant
    # Global wavelet spectrum & significance levels: 
    #same but for global calclation we are
    # doing a time-averaged significance test - hence we pass `sigtest=1`
    dof = num - scale  # the -scale corrects for padding at edges
    global_signif=wave_signif(variance,dt=dt,scale=scale,sigtest=1,dof=dof,mother=mother)
    glsig95 = global_signif[:, np.newaxis].dot(np.ones(num)[np.newaxis, :]) 
        # expand signif --> (J+1)x(N) array
    normsig= power/max(map(max,power))
    normsig= 1- normsig
    return time_array, period, scale, coi, power, global_ws, sig95, glsig95, normsig, global_signif



##############################################################


