import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib.dates import DateFormatter


def spectrum_plotter(SAVEPATH, title, t, data_shift, env, data_shift_norm, power, period, sig95, coi, global_signif, global_ws, errorbar_packed):

    solarmax, solarmin, period_max_power_overall, min_sig95_maxpow, max_sig95_maxpow, index_max_power \
    ,period_max_power_solarmax, min_sig95_solarmax,max_sig95_solarmax, index_solarmax = errorbar_packed
    datatype=r"Frequency shift ($\mu$Hz)" 
    cl_sm='dodgerblue'
    cl_maxp='red'

    fig = pylab.figure(figsize = (11,8))
    ax = pylab.axes([0.069, 0.7, 0.68, 0.25]) #top
    dx = pylab.axes([0.069, 0.45, 0.68, 0.25]) #middle
    bx = pylab.axes([0.069, 0.1, 0.68, 0.35]) #bottom
    cx = pylab.axes([0.755, 0.1, 0.19, 0.35], sharey=bx) #bottom right

    ##### First plot the input timeseries ####

    fig.suptitle(title, fontsize=16)
    ax.plot(t, data_shift, 'k', label='Data')
    ax.plot(t, env, 'r', label='Envelope')
    ax.set_ylabel(datatype)
    ax.grid(which='major', axis='x')
    ax.legend(loc='upper right')

    dx.plot(t, data_shift_norm, 'k')
    dx.set_ylabel('Detrended \n '+datatype)
    dx.grid(which='major', axis='x')

    #ax.set_xlim(tnew2[0], tnew2[-1])
    ax.tick_params(axis="x", which="both", labelbottom=False)
    dx.tick_params(axis="x", which="both", labelbottom=False)

    #### Plot the wavelet power spectrum #####
    levels = np.linspace(0, np.max((power)), 20)
    bx.set_yscale("linear")
    cs = bx.contourf(t, period, power, cmap="inferno", levels=levels, extend='both')
    contours = bx.contour(t, period, sig95, [1.0], extend='both',colors='white', linewidths=2.)
    bx.fill_between(t, coi, np.full((len(t)), 5000), linewidth=0, facecolor='white', edgecolor="white", hatch="X", alpha=0.25)
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

    ############
    #ERRORBARS
    ############

    #Overall_power
    bx.plot(t[index_max_power],period_max_power_overall,color=cl_maxp,marker='*',markersize=10)
    bx.plot([t[index_max_power], t[index_max_power]], [min_sig95_maxpow,max_sig95_maxpow], color=cl_maxp, linewidth=1.5)
 
    #SolarMax
    #bx.plot([t[index_solarmax], t[index_solarmax]], [min_sig95_solarmax,max_sig95_solarmax], color=cl_sm, linewidth=1.5)
    #bx.plot(t[index_solarmax], period_max_power_solarmax ,color=cl_sm,marker='*',markersize=10)

    #bx.legend([plt.Line2D((0,1),(0,0), color=cl_maxp, marker='*', linestyle=''),plt.Line2D((0,1),(0,0), color=cl_sm, marker='*', linestyle='')],['Maximal power','Maximal power at Solar Max'],loc='upper right')
    bx.legend([plt.Line2D((0,1),(0,0), color=cl_maxp, marker='*', linestyle='')],['Maximal power'],loc='upper right')

    cx.set_xscale('log')
    cx.set_yscale('log')
    cx.tick_params(axis='y', which='major', labelleft=False)

    plt.savefig(SAVEPATH+'/'+title+'.png', dpi=200)
    plt.close()

