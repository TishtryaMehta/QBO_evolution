import numpy as np
import matplotlib.pyplot as plt
import datetime
from datetime import date



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
    if overwrite< 20:
        p = contour_data.collections[0].get_paths()[overwrite]
        v = p.vertices
        t_contour = v[:,0]
        per_contour = v[:,1]
    elif overwrite >=20:
        length_contour=[]
        for n in np.arange(0, len(contour_data.collections[0].get_paths())):
            #print('n is', n)
            p = contour_data.collections[0].get_paths()[n]
            v = p.vertices
            t_contour = v[:,0]
            per_contour = v[:,1]
            length_contour.append(len(t_contour))
            
        # This selects the longest contour; this isn't always the correct contour
        # You may manually need to select the correct contour
        # by changing the value of max_length_contour to any number
        # in n until it fulfils your requirements 
        max_length_contour=length_contour.index(max(length_contour))
        
        # UNCOMMENT TO VISUALISE AND CHOOSE CONTOUR
        #print('Lengths of the contours: ', length_contour, )
        #print('The index of the largest contour is: ', max_length_contour)
      
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
        #print('There are no upper and lower periods for this index')
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
        solarmax_period_min=periods_contours_index_array[temp_diff_index[0][0]]
    else: 
        solarmax_period_min=min(period_data)
    ### MAXIMUM PERIOD ###
    if len(diff_periods[diff_periods>0]) != 0:
        temp_diff_index=np.where(diff_periods==min(diff_periods[diff_periods>0]))
        solarmax_period_max=periods_contours_index_array[temp_diff_index[0][0]]
    else:
        solarmax_period_max=max(period_data)
    
    return period_max_power, solarmax_period_min, solarmax_period_max, t_contour, per_contour


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
# PERIOD OVERALL MAX POWER MIN AND MAX 95% SIG
# Finds the minimum and maximum period which corresponds 
# to 95% significance contour at the index of the maximum power 
#""""""""""""""""""""""""""""""""""""""""""""""""""""

def minmax_maxpowerfinder(contours, max_power_index, period, power,time,overwrite):
    epoch_date=date(1970,1,1)
    max_power_index_date=date(time[max_power_index].year,time[max_power_index].month,time[max_power_index].day)
    delta_days_temp=max_power_index_date-epoch_date
    delta_days_temp_days=delta_days_temp.days
    temp, minimum_sig95,maximum_sig95, temp, temp=minmaxperiod(contours,max_power_index, period, power, delta_days_temp_days,overwrite)
    return minimum_sig95,maximum_sig95

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

  
    
def errorbar_finder(t, data_raw, power, period, sig95, seg1_index_min, seg1_index_max):
    
    #obtain overall max power and period 
    max_power_period, index_max_power=period_overall_max_power(power, period,seg1_index_min,seg1_index_max)
    fig, ax = plt.subplots()
    contours = ax.contour(t, period, sig95, [1.0])
    total_contours=len(contours.collections[0].get_paths())

    plt.close()

    #Obtain Solarmax, Solarmin and index_solarmax
    solarmax, solarmin, delta_days, index_solarmax=solarmax_finder(t,data_raw, seg1_index_min, seg1_index_max)
    #Obtain solarmax_period, solarmax_period_min,solarmax_period_max,

    solarmax_period, solarmax_period_min,solarmax_period_max, t_contour, per_contour=minmaxperiod(contours,index_solarmax, period, power, delta_days,0)

    #Obtain max_power_min,max_power_max
    minimum_lower_error_diff=[]
    maximum_upper_error_diff=[]
    for contour_chooser in np.arange(0,total_contours):
        max_power_min,max_power_max=minmax_maxpowerfinder(contours, index_max_power, period, power,t,contour_chooser)
        minimum_lower_error_diff.append(max_power_period-max_power_min)
        maximum_upper_error_diff.append(max_power_max-max_power_period)
        #print(max_power_min, contour_chooser)

    max_power_max=max_power_period + min(maximum_upper_error_diff)
    max_power_min= max_power_period- min(minimum_lower_error_diff)

    return solarmax, solarmin, max_power_period, max_power_min, max_power_max, index_max_power, solarmax_period, solarmax_period_min, solarmax_period_max, index_solarmax
