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
# WRITTEN BY T MEHTA, 2023
# CONTACT TISHTRYA.MEHTA@GMAIL.COM IF YOU HAVE ANY PROBLEMS :)

# THIS CODE USES THE TORRENCE & COMPO WAVELET CODE, IN `waveletFunctions.py`
# AND OUTPUTS THE SIGNIFICANT PERIODS AND ERRORS IN PERIODICITY FOR THE
# QBO USING HELIOSESIMIC AND RADIO DATA
# CODE IS WRITTEN FOR PYTHON 3.6.9

# INPUTS:
# Source == "MDIHMI" or "GONG" or "F107"
# inputcycles == "23","24" or "2324"
# inputfreq == "high", "low" or "all"


def QBO_cwt_finder(datasource, inputcycles, inputfreq):


    ##############################################################
    # Read in data and sort according to datatype and duration
    ##############################################################

    # General path containing directories of data
    path="/home/space/phrsnz/Desktop/PhD/QBO_evolution/Data_processed/"
    filetype=".dat"
    #Column where the data is saved, default is zeroth column
    index=0
    #Type of data, default is frequency shift
    datatype=r"Frequency shift ($\mu$Hz)" 
    cl_sm='dodgerblue'
    cl_maxp='red'

    
    # MDIHMI DATA
    if str(datasource)==str('MDIHMI'):
        dt=72 #days
        if str(inputcycles)==str('23') or str(inputcycles)==str('2324'):
            t1=datetime.datetime(1996,7,12)
        elif str(inputcycles)==str('24'):
            t1=datetime.datetime(2008,12,12)

    # GONG DATA
    elif str(datasource)==str('GONG'):
        dt=36 #days
        if str(inputcycles)==str('23') or str(inputcycles)==str('2324'):
            t1=datetime.datetime(1996,9,3) #YEAR MONTH DAY 
        elif str(inputcycles)==str('24'):
            t1=datetime.datetime(2008,12,29)

    # F10.7 DATA
    elif str(datasource)==str('F107'):
        filetype=".txt"
        index=3 #column where the data is saved
        dt=1 #days
        datatype="10.7 cm Radio Flux (sfu)"
        if str(inputcycles)==str('23') or str(inputcycles)==str('2324'):
            t1=datetime.datetime(1996,8,1) #YEAR MONTH DAY 
        elif str(inputcycles)==str('24'):
            t1=datetime.datetime(2008,12,1)

    datapath=str(datasource)+"_cycle_"+str(inputcycles)+"/"+str(datasource)+"_"+str(inputfreq)+"_"+str(inputcycles)+filetype
    savepath="/home/space/phrsnz/Desktop/PhD/QBO_evolution/CWT_analysis/"
    data_full=np.genfromtxt(path+datapath)
    data=[item[index] for item in data_full]
    data_raw=np.copy(data)

    ##############################################################
    # DEFINE FUNCTIONS
    ##############################################################

    #""""""""""""""""""""""""""""""""""""""""""""""""""""
    # MINMAXPERIOD
    # Outputs min and max periods correpsonding to 
    # the confidence level for a given index (e.g. Solar Max)
    #""""""""""""""""""""""""""""""""""""""""""""""""""""
    def minmaxperiod(contour_data,index_data, period_data, power_data,delta_days, overwrite):
        #Obtain period corresponding to maximum power in given index
        period_max_power= period_data[np.where(power_data[:,index_data]==max(power_data[:,index_data]))][0]
        
        # We can easily obtain the periods that are above 1 in sig95 by printing a masked
        # array of the sig95 array with values >1 and printing the corresponding periods
        # However I want to find the period correspnding to the contours of the confidence level so I obtain the 
        # paths of the contours: 

        ### SELECT CORRECT CONTOUR PATH ###
        #In some cases, the contour data is not a list of 1 paths, but a list of 2 or more paths.
        #in this case, we take the longest contour (with the largest perimeter)
        # for all the cases do p =  contour_data.collections[0].get_paths()[n] from n=0 to len(get_paths)
        if overwrite ==1: #Overwrite to select specific paths
            p = contour_data.collections[0].get_paths()[3] #For all freq G2324 we want the 3rd contour
            v = p.vertices
            t_contour = v[:,0]
            per_contour = v[:,1]
        elif overwrite ==2: #Overwrite to select specific paths
            p = contour_data.collections[0].get_paths()[4] # THIS WAS 4 change BACK for not MDIHMI 2324 HIGH
            v = p.vertices
            t_contour = v[:,0]
            per_contour = v[:,1]
        else:
            length_contour=[]
            for n in np.arange(0, len(contour_data.collections[0].get_paths())):
                #print('n is', n)
                p = contour_data.collections[0].get_paths()[n]
                v = p.vertices
                t_contour = v[:,0]
                per_contour = v[:,1]
                length_contour.append(len(t_contour))

            max_length_contour=length_contour.index(max(length_contour))
            #refind the data for the largest contour
            p = contour_data.collections[0].get_paths()[max_length_contour]
            v = p.vertices    
            t_contour = v[:,0]
            per_contour = v[:,1]

        ### SELECT UPPER AND LOWER PERIODS ###
        #Find indexes where t_contour (the time coordinate of the contours) is equal to the time coordinate of solar max
        #This should give 2 (possibly more) indexes corresponding to the upper and lower periods of the confidence level
        #The time coordniates for t and t_contour are not the same; t is in datetime format and t_contour is in days since
        # 1st Jan 1970 (Unix time). So we find Delta days between solar max and 1st Jan 1970 and find the index of t_contour
        # where the time coordinate is equal to this delta day
        min_max_indexes=np.where(np.round(t_contour) == delta_days)
        #print(t_contour, 't_contour')
        #print(delta_days, 'delta_days')
        #print(min_max_indexes, 'min_max_indexes')
        if len(min_max_indexes[0])==0:
            print('No contour found')
            return period_max_power,min(period_data),max(period_data),t_contour, per_contour
        #find the periods corresponding to all the values in min_max_indexes
        periods_contours_index_array=[]
        for b in np.arange(0, len(min_max_indexes[0])):
            periods_contours_index=per_contour[min_max_indexes[0][b]]
            periods_contours_index_array.append(periods_contours_index)
        #of the values in periods_contours_index_array find the value closest to the period corresponding to the maximum
        #power that is greater than the period corresponding to the maximum power calculate the difference between 
        #periods_contours_index_array and period_max_power
        diff_periods=np.array(periods_contours_index_array)-period_max_power
        #find the index of the value in diff_periods that is the smallest positive number
        ### MINIMUM PERIOD ###
        if len(diff_periods[diff_periods<0]) != 0:
            temp_diff_index=np.where(diff_periods==max(diff_periods[diff_periods<0]))
            min_sig95_solarmax=periods_contours_index_array[temp_diff_index[0][0]]
        else: 
            min_sig95_solarmax=min(period_data)
        ### MAXIMUM PERIOD ###
        if len(diff_periods[diff_periods>0]) != 0:
            temp_diff_index=np.where(diff_periods==min(diff_periods[diff_periods>0]))
            max_sig95_solarmax=periods_contours_index_array[temp_diff_index[0][0]]
        else:
            max_sig95_solarmax=max(period_data)
        
        return period_max_power, min_sig95_solarmax, max_sig95_solarmax, t_contour, per_contour

    #""""""""""""""""""""""""""""""""""""""""""""""""""""
    # PERIOD OVERALL MAX POWER
    # Finds the period which corresponds 
    # to the overall maximum power in the data 
    # between indexmin and indexmax
    #""""""""""""""""""""""""""""""""""""""""""""""""""""

    def period_overall_max_power(power_data, period_data, indexmin,indexmax):
        #Find the index of the maximum power
        #print('shape power',np.shape(power_data[:,70:-1]))
        #print(power_data[:,indexmin:indexmax] )
        index_max_power=np.where(power_data[:,indexmin:indexmax] == max(map(max,power_data[:,indexmin:indexmax])))[1][0]
        index_max_power=np.copy(index_max_power)+indexmin
        #print('index_max_power',index_max_power)
        #print('time:', t[index_max_power])
        #print('power:', power_data[:,index_max_power])
        #print('period:', period_data[np.where(power_data[:,index_max_power]==max(power_data[:,index_max_power]))][0])
        #Find the period corresponding to the maximum power
        period_max_power= period_data[np.where(power_data[:,index_max_power]==max(power_data[:,index_max_power]))][0]
        return period_max_power, index_max_power


    #""""""""""""""""""""""""""""""""""""""""""""""""""""
    # SOLARMAX FINDER
    # Finds the dates of solar max and solar min defined 
    # as the min and max of the input data within the 
    # indexes prescribed
    #""""""""""""""""""""""""""""""""""""""""""""""""""""

    def solarmax_finder(time_data,input_data, lower_index, upper_index):
        ###POWER SOLARMAX AND SOLARMIN
        #Define epoch_date as 1st Jan 1970
        epoch_date=date(1970,1,1)
        #define the max_amp as the maximum amplitude for data_raw
        input_data=np.copy(input_data[lower_index:upper_index])
        max_amp=max(input_data)
        #find the index of data_raw that corresponds to its maximual value
        index_max_amp=np.where(input_data == max(input_data))[0]
        index_max_amp=np.copy(index_max_amp[0])
        #define index_min_amp as the minimum amplitude for data_raw
        index_min_amp=np.where(input_data == min(input_data))[0]
        index_min_amp=np.copy(index_min_amp[0])
        #define solarmin as the time corresponding to index_min_amp
        solarmin=time_data[index_min_amp+lower_index]
        #define solarmax as the time corrsponding to index_max_amp
        solarmax=time_data[index_max_amp+lower_index]
        #Convert Solar max to a date object
        solarmax_date=date(solarmax.year,solarmax.month,solarmax.day)
        #Calculate the difference between solar max and epoch date in days
        delta=solarmax_date-epoch_date
        #Define delta_days to be carried in minmaxperiod function
        delta_days=delta.days
        return solarmax, solarmin, delta_days, index_max_amp+lower_index

    ##############################################################

    
    #Make sure the length of the data is even (needed for the below wavelet analysis)   
    if (len(data_raw)% 2) != 0:
        data_raw=np.copy(data_raw[:-1])
  
    # Define the envelope of the data
    if str(inputcycles)!=str('2324'):
        env=savgol_filter(data_raw,len(data_raw)-31,5) #-31, 5    #savgol
    else:
        env=savgol_filter(data_raw,len(data_raw)-40,7) #FOR GONG 2324
    #env=ndimage.median_filter(data_raw, size=int(len(data_raw)/5), mode='nearest')#) #/5 good for f10.7  #medin filter
    #env= spline(data_raw, order=2, dspline=1000) #spline interpolation
    #Subtract the envelope from the data
    signal2=data_raw-env

    #Set up the time axis
    t=[t1 + datetime.timedelta(days =dt*i) for i in range(len(signal2))]
    
    num=len(signal2)
    time=np.arange(num)*dt
    mother = 'MORLET' 
    # here we'll use a Morlet mother wavelet, 
    #this can be changed see notes in waveletFunctions.py
    variance = np.std(signal2, ddof=1) ** 2 # variance of the signal. 

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
    wave, period, scale, coi = wavelet(signal2, dt, mother="MORLET")
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
    lag1=0.72#2.1
    signif = wave_signif(variance, dt=dt, sigtest=0, scale=scale,lag1=lag1, mother=mother) 
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


    """------------------------------
    Plot the results!
    -------------------------------"""

    fig = pylab.figure(figsize = (11,8))
    ax = pylab.axes([0.069, 0.7, 0.68, 0.25]) #top
    dx = pylab.axes([0.069, 0.45, 0.68, 0.25]) #middle
    bx = pylab.axes([0.069, 0.1, 0.68, 0.35]) #bottom
    cx = pylab.axes([0.755, 0.1, 0.19, 0.35], sharey=bx) #bottom right

    ##### First plot the input timeseries ####

    fig.suptitle(str(datasource)+ ' Cycle(s) '+str(inputcycles)+' over '+str(inputfreq)+' frequencies')
    ax.plot(t, data_raw, 'k', label='Data')
    ax.plot(t, env, 'r', label='Envelope')
    ax.set_ylabel(datatype)
    ax.grid(which='major', axis='x')
    ax.legend(loc='upper right')

    dx.plot(t, signal2, 'k')
    dx.set_ylabel('Detrended \n '+datatype)
    dx.grid(which='major', axis='x')

    #ax.set_xlim(tnew2[0], tnew2[-1])
    ax.tick_params(axis="x", which="both", labelbottom=False)
    dx.tick_params(axis="x", which="both", labelbottom=False)

    #### Plot the wavelet power spectrum #####
    levels = np.linspace(0, np.max((power)), 20) #17
    # you can change these levels to highlight the power
    bx.set_yscale("linear")
    cs = bx.contourf(t, period, power, cmap="inferno", levels=levels, extend='both')
    # plot the power
    #print(0.95*max(d1))
    contours = bx.contour(t, period, sig95, [1.0], extend='both',colors='white', linewidths=2.)


    # plot the significance level contours
    #bx.plot(t, coi, 'w--', alpha=0.25) # plot the cone of influence

    bx.fill_between(t, coi, np.full((len(t)), 5000), linewidth=0, facecolor='white', edgecolor="white", hatch="X", alpha=0.25)
    #bx.legend(loc='upper right')
    bx.set_ylabel('Period (Days)')
    bx.set_xlabel('Years')
    bx.invert_yaxis()
    bx.set_ylim(100,4000)
    bx.set_xlim(t[0], t[-1])
    bx.xaxis.set_major_formatter(DateFormatter('%Y'))
    bx.fmt_xdata=DateFormatter('Y:%')
    fig.autofmt_xdate() 
    #### Plot the global wavelet spectrum #####
    cx.plot(global_signif, (period), color = 'r', label = 'Sig 95%', lw = 1.5)
    cx.plot(global_ws, period, 'k-', linewidth=1.5)
    #cx.legendloc = ('upper right')
    cx.set_title('Global Wavelet')
    cx.set_xlabel(r'Power')


    # If only one cycle, find overall max power. If two cycles perform operations twice:
    if inputcycles != '2324':

        ###POWER OVERALL
        seg1_index_min,seg1_index_max = 0,-1
        #Obtain the period corresponding to the overall maximum power
        period_max_power_overall, index_max_power=period_overall_max_power(power, period,seg1_index_min,seg1_index_max)
        #plot the period corresponding to the maximum power on the power spectrum with a red star
        bx.plot(t[index_max_power],period_max_power_overall,color=cl_maxp,marker='*',markersize=10)

        #Obtain Solarmax, solarmin and delta_days
        solarmax, solarmin, delta_days, index_solarmax=solarmax_finder(t,data_raw, seg1_index_min, seg1_index_max)

        #Obtain min max periods for solar max
        period_max_power, min_sig95_solarmax,max_sig95_solarmax, t_contour, per_contour=minmaxperiod(contours,index_solarmax, period, power, delta_days,0)

        # Plot a line between min_sig95_solarmax and max_sig95_solarmax in pink
        bx.plot([t[index_solarmax], t[index_solarmax]], [min_sig95_solarmax,max_sig95_solarmax], color=cl_sm, linewidth=1.5)
        #plot the maximal power correponding to the maximum power at time t[index_max_amp] on the power spectrum with a pink star
        bx.plot(t[index_solarmax], period_max_power ,color=cl_sm,marker='*',markersize=10)


    else:
        seg1_index_min,seg1_index_max = 0,round(len(data_raw)/2)
        seg2_index_min,seg2_index_max = round(len(data_raw)/2),len(data_raw)
        solarmax_1, solarmin_1, delta_days_1, index_solarmax1 =solarmax_finder(t,data_raw, seg1_index_min, seg1_index_max)

        solarmax_2, solarmin_2, delta_days_2, index_solarmax2=solarmax_finder(t,data_raw, seg2_index_min, seg2_index_max)
        period_max_power_seg1, index_max_power_seg1=period_overall_max_power(power, period,seg1_index_min,seg1_index_max)
        #print(period_max_power_seg1, index_max_power_seg1)
        period_max_power_seg2, index_max_power_seg2=period_overall_max_power(power, period,seg2_index_min,seg2_index_max)
        #Obtain min max periods for solar max
        period_solarmax_1, min_sig95_solarmax_1,max_sig95_solarmax_1, t_contour, per_contour=minmaxperiod(contours,index_max_power_seg1, period, power, delta_days_1,0)
        bx.plot([t[index_solarmax1], t[index_solarmax1]], [min_sig95_solarmax_1,max_sig95_solarmax_1], color=cl_sm, linewidth=1.5)
        bx.plot(t[index_solarmax1],period_solarmax_1,color=cl_sm,marker='*',markersize=10)
        #set overwrite_val to equal 1 so force the min and max periods to be calculated for the second contour
        # he below is finickity and needs to be changed for different datasources as each wavelet spectra 
        #produces a different number of contours. I have manually tuned which contour is selected
        # for the GOES and MDI HMI data. Overwrite_val tells the minmaxperiod function which contour to use. 

        if str(datasource) == 'GONG':
                overwrite_val=1
        elif str(datasource) == 'MDIHMI' and str(inputfreq)!='high':
                overwrite_val=2
        elif str(datasource) == 'MDIHMI' and str(inputfreq) =='high':
                overwrite_val=1

        period_solarmax_2, min_sig95_solarmax_2,max_sig95_solarmax_2, t_contour, per_contour=minmaxperiod(contours,index_max_power_seg2, period, power, delta_days_2,overwrite_val)
        bx.plot([t[index_solarmax2], t[index_solarmax2]], [min_sig95_solarmax_2,max_sig95_solarmax_2], color=cl_sm, linewidth=1.5)
        bx.plot(t[index_solarmax2],period_solarmax_2,color=cl_sm,marker='*',markersize=10)

        #plot the period corresponding to the maximum power on the power spectrum with a red star
        bx.plot(t[index_max_power_seg1],period_max_power_seg1,color=cl_maxp,marker='*',markersize=10)
        bx.plot(t[index_max_power_seg2],period_max_power_seg2,color=cl_maxp,marker='*',markersize=10)

        #NEED TO WRITE CODE TORECOGNISE EACH SOLARE MAX IN THE DATA HERE AND BASICALLY REPEAT THE PROCEDURES
        #FOR TWO SOLAR MAXES :) 



    #Add a legend to bx with red star symbol  = true maximum power and pink star symbol = maximum power at solarmax
    bx.legend([plt.Line2D((0,1),(0,0), color=cl_maxp, marker='*', linestyle=''),plt.Line2D((0,1),(0,0), color=cl_sm, marker='*', linestyle='')],['Maximal power','Maximal power at Solar Max'],loc='upper right')
    ###overplot linfitcomp
    cx.set_xscale('log')
    cx.set_yscale('log')
    #cx.set_ylim(0,1000)
    #cb = plt.colorbar(cs, orientation = 'vertical') 
    # colorbar for the power plot
    #cb.set_label('Wavelet power')
    cx.tick_params(axis='y', which='major', labelleft=False)
    
    plt.savefig(savepath+str(datasource)+"_"+str(inputcycles)+"_"+str(inputfreq)+".png", dpi=200)

    #""""""""""""""""""""""""""""""""
    # Print Output
    #""""""""""""""""""""""""""""""""
    
    if inputcycles != '2324':
        print("******************************************************************************")
        print("............................................")
        print("You have run the QBO finder code for ",datasource," data, over ",inputfreq, "frequencies and ",inputcycles, " solar cycle(s).")
        print("Your output is saved here:", savepath, 'as:', str(datasource),"_",str(inputcycles),"_",str(inputfreq),".png")
        print("............................................")
        print('Solarmax: ',solarmax, ', Solarmin: ',solarmin)
        print("Period with max power at Solar max: %.2f d"%period_max_power)
        print('Min Period at CI at Solar max: %.2f d'%min_sig95_solarmax)
        print('Max Period at CI at Solar max: %.2f d'%max_sig95_solarmax)
        print('Periodicity with maximum power between time %s -- %s : %.2f d'%(t[seg1_index_min],t[seg1_index_max],period_max_power_overall))
        print("******************************************************************************")

        csvname='Single_cycle_QBO_output.csv'
        csvpath= "/home/space/phrsnz/Desktop/PhD/QBO_evolution/CWT_analysis/Output/"
    
        header="#Instrument, #Input Cycle, # Input Freq, # Solar Max, # Solar Min, # Period Max Power, # Min Period CI, # Max Period CI, # Period Max Power Overall"
        data_list = [datasource, inputcycles, inputfreq, solarmax, solarmin, period_max_power, min_sig95_solarmax, max_sig95_solarmax, period_max_power_overall]

        file_exists = os.path.isfile(csvpath+csvname)

        with open (csvpath+csvname, 'a') as f_object:
            if not file_exists:
                writer_object = writer(f_object)
                writer_object.writerow([header])
   
            writer_object = writer(f_object)
            writer_object.writerow(data_list)
            f_object.close()


    else:
        print("******************************************************************************")
        print("............................................")
        print("You have run the QBO finder code for ",datasource," data, over ",inputfreq, "frequencies and ",inputcycles, " solar cycle(s).")
        print("Your output is saved here:", savepath, 'as:', str(datasource),"_",str(inputcycles),"_",str(inputfreq),".png")
        print("............................................")
        print('Solarmax (Cycle 23): ',solarmax_1, ', Solarmin (Cycle 23): ',solarmin_1)
        print('Solarmax (Cycle 24): ',solarmax_2, ', Solarmin (Cycle 23): ',solarmin_2)
        print('#####CYCLE 23#####')
        print("Period with max power at Solar max (23): %.2f d"%period_max_power_seg1)
        print('Min Period at CI at Solar max: %.2f d'%min_sig95_solarmax_1)
        print('Max Period at CI at Solar max: %.2f d'%max_sig95_solarmax_1)
        print('Periodicity with maximum power between time %s -- %s : %.2f d'%(t[seg1_index_min],t[seg1_index_max],period_max_power_seg1))
        print('#####CYCLE 24#####')
        print("Period with max power at Solar max (24): %.2f d"%period_max_power_seg2)
        print('Min Period at CI at Solar max: %.2f d'%min_sig95_solarmax_2)
        print('Max Period at CI at Solar max: %.2f d'%max_sig95_solarmax_2)
        print('Periodicity with maximum power between time %s -- %s : %.2f d'%(t[seg2_index_min],t[seg2_index_max-1],period_max_power_seg2))
        print("******************************************************************************")

        csvname='Two_cycle_QBO_output.csv'
        csvpath= "/home/space/phrsnz/Desktop/PhD/QBO_evolution/CWT_analysis/Output/"
    
        header="#Instrument, #InputCycle, # InputFreq, #T_SolarMax23, #T_SolarMin23, #PeriodMaxPower23, #MinPeriodCI23, # MaxPeriodCI23, # PeriodMaxPowerOverall23, #T_SolarMax23, #T_SolarMin24, #PeriodMaxPower24, #MinPeriodCI24, # MaxPeriodCI24, # PeriodMaxPowerOverall24"
        data_list = [datasource, inputcycles, inputfreq, solarmax_1, solarmin_1, period_max_power_seg1, min_sig95_solarmax_1, max_sig95_solarmax_1, period_max_power_seg1,solarmax_2, solarmin_2, period_max_power_seg2, min_sig95_solarmax_2, max_sig95_solarmax_2, period_max_power_seg2]


        file_exists = os.path.isfile(csvpath+csvname)

        with open (csvpath+csvname, 'a') as f_object:
            if not file_exists:
                writer_object = writer(f_object)
                writer_object.writerow([header])
   
            writer_object = writer(f_object)
            writer_object.writerow(data_list)
            f_object.close()



#*************************************************************************************
#Run code

for a in ['GONG','MDIHMI']:
    for b in ['23', '24', '2324']:
        for c in ['low','high','all']:
            x=QBO_cwt_finder(a,b,c)

