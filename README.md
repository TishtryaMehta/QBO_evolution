# QBO_evolution
Data and codes for investigating the evolution of the instantaneous period of the solar QBO in helioseismic data. 

The technique can be applied to any quasi-oscillatory data to extract statistically significant oscillations, aobut the 98% Confidence Level. 

A full user guide to this code, including examples of outputs, can be found here: https://docs.google.com/document/d/1KfyH5kAgkOeTsG2qZ-GQ2lZDPcurB6GKnvs8vY1TYgo/edit?usp=sharing

This code was written for Python 3.9.6, July 2023.
Any questions can be sent to: tishtrya.mehta@gmail.com

________________________________________________________________________________________________________________________________

Below is an abridged excerpt of the user guide, describing how to use wavelet_analysis.py
________________________________________________________________________________________________________________________________


Section 1. Performing Wavelet analysis on user inputted data, using wavelet_analysis.py
Wavelet_analysis.py allows a user to read in data in csv or txt format and extract the significant periodicity from the normalised data’s continuous wavelet transform (CWT) and Global Wavelet Transform (GWT) with error bars corresponding to  98% confidence levels. 

The code works as follows:

The data is read in using  wavelet_data_reader , and reads the column indicated by variable index,reading in the filename of the data, where the data has been saved in DATA_PATH. 
The data is preprocessed by removing the nan values from the dataif necessary removing the final data point to ensure the data is even
Normalising the data by finding the best fitting savgol envelope for the data, and dividing the data by it. The savgol envelope can be manipulated by changing the variables savgol_wl and savgol_pol. This should be done manually, by eye, until a good fit is found. A good fit should result in the normalised data having a roughly constant amplitude (Shown in the second panel of the produced figure).
This process uses cadence dt, and start time starttime.
Note: We normalise the data because we may want to only examine the overall periodicities in the data and ignore the effects of amplitude modulation. Wavelet analysis is sensitive to the amplitude of the data by design and assessing non-normalised data may mean you miss observing periodicities at low amplitudes/ your signal may be dominated by large scale oscillations. However if you wish to examine amplitude modulation you can force your savgol envelope to be constant ([1,1,1,..]) which will remove the normalisation effect. 
The normalised data is assessed by wavelet_analysis_finder , producing a CWT, GWT and associated coi, etc, 
We find the period corresponding to the maximum power in the CWT. We also obtain the corresponding errors on this period. We do this by assessing the 98% confidence contour, and finding the minimum and maximum periods that correspond to the 98% contour at the time of the maximum power period. We only search over the range specified by the user inputted variables seg1_index_min, seg1_index_max. 
The above is then plotted using  spectrum_plotter 
The above is then saved as a figure, with title as a user defined variable 

Alongside downloading the input data, the user will have to change the following variables:

#    User chosen variables
      ----------
      filename : string
          Name of the file containing the data (csv or txt format)
      index : int
          The index of where the data is saved in the input file
      savgol_wl : int
          The window length for the savgol filter
      savgol_pol : int
          The polynomial order for the savgol filter
      starttime : datetime.datetime
          The start time of the data (YYYY, M, D)
      dt : int
          The time step of the data in days
      title : string
          The title of the output plot and csv file
      seg1_index_min : int                                            ****
          The minimum index of the data to be examined
      seg1_index_max : int
          The maximum index of the data to be examined
      DATA_PATH : string
          The path to where the input data is saved
      SAVE_PATH : string
          Where the output plot and csv file are saved

**** To examine the full duration of the data seg1_index_min should =0 and seg1_index_max should = -1. 

E.g.:
filename="v1f_1000-2000.csv"
savgol_wl=31
savgol_pol=5
starttime=datetime.datetime(1996,9,3)
dt=36
title='All_modes_1000-2000uHz'
seg1_index_min,seg1_index_max = 0,-1
index=0
DATA_PATH =  os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "CWT_analysis_fmodes", "Output")
sAVE_PATH =  os.path.join("/", "home", "space", "phrsnz", "Desktop", "Academic", "Programs","PhD", "QBO_evolution", "CWT_analysis_fmodes", "Output")
#---------------------------

After implementing these variables to your preference, the code can be run. It relies on

V1f_data_sorting.py
Wavelet_spectrum_maker.py
Wavelet_errorbars.py
Wavelet_spectrum_plotter.py
Wavelet_data_saver.py
Constants.py

which work in the order above to read in the data and preprocess the data, perform wavelet analysis on the preprocessed data, find the corresponding errors, create a plot/ figure of the analysis, and save the spectrum and the results into a csv file. 

CSV:
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
e.g.
"T_[MaxAmp]","T_[MinAmp]","max_pow_period","min_pow","max_pow","t_max_pow_period","period_solar_max","min_pow_solar_max","max_pow_solar_max"
2001-09-13 00:00:00, 2007-07-08 00:00:00, 637.74, 269.88, 877.36, 2001-04-22 00:00:00,623.17,74.38,4441.48

