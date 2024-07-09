import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import iqr
import pandas as pd
from scipy import signal, stats
from matplotlib import mlab
import os
import sys
func_path = os.path.join(os.path.dirname(__file__), '../required_functions')
sys.path.append(func_path)
    
##path to pre processed surface data
path_to_surface_data='../Pre_processed_data'
## Load surface data
upper=np.load(path_to_surface_data+'/upper.npy')
lower=np.load(path_to_surface_data+'/lower.npy')
# ap=upper-lower

## quick way to find H and Sp using filtering
def FFT_H_int_all_traces_one_plot(trace_in):
    trace = signal.detrend(trace_in, axis=0)
    narray = trace.shape[0]
    trace = trace.flatten()
    Pxx, _ = mlab.psd(trace, NFFT=narray-1)
    kvals = np.arange(1, narray // 2 + 1)
    return kvals, Pxx

def calculate_scaling_PSD(surface):
    x_length, y_length = surface.shape

    def process_traces(traces, length):
        lower_kvals = []
        lower_Pxx = []
        for i in range(length):
            t = traces[i, :].reshape(length, 1)
            x, y = FFT_H_int_all_traces_one_plot(t)
            lower_kvals.append(x)
            lower_Pxx.append(y)
        return lower_kvals, lower_Pxx

    lower_kvals_x, lower_Pxx_x = process_traces(surface, x_length)
    lower_kvals_y, lower_Pxx_y = process_traces(surface.T, y_length)

    comb_px = np.vstack((lower_Pxx_x, lower_Pxx_y))
    Px = np.mean(comb_px, axis=0)
    kvals = lower_kvals_x[0]

    vals_to_remove = 1
    init_rem = 2

    Px_mean_no_first = Px[init_rem:-vals_to_remove]
    kvals_mean_no_last = kvals[:-(vals_to_remove + init_rem)]

    log_kvals = np.log10(kvals_mean_no_last)
    log_Px = np.log10(Px_mean_no_first)
    res = stats.linregress(log_kvals, log_Px)

    slope = res.slope
    c = 10**res.intercept
    H_out = (slope + 1) / -2
    K = kvals_mean_no_last
    N = (len(kvals_mean_no_last) * 2)
    a = ((2 * np.sqrt(2)) / N) * np.sqrt(c)
    N = (N / 2) - 1
    b = np.sqrt(np.sum((K**-(H_out + 0.5)) * np.sin(np.pi * (K / N)))**2)
    std_h = a * b
        
    return H_out, std_h

## lower
H_lower, Sp_lower=calculate_scaling_PSD(lower)
print (H_lower, Sp_lower)
np.save('calc_H_PSD_unfiltered/H_lower.npy',H_lower)
np.save('calc_H_PSD_unfiltered/Sp_lower.npy',Sp_lower)

## upper
H_upper, Sp_upper=calculate_scaling_PSD(upper)
print (H_upper, Sp_upper)
np.save('calc_H_PSD_unfiltered/H_upper.npy',H_upper)
np.save('calc_H_PSD_unfiltered/Sp_upper.npy',Sp_upper)




def plot_unfiltered_PSD(surface,fig_lab, output):
    x_length, y_length = surface.shape

    def process_traces(traces, length):
        lower_kvals = []
        lower_Pxx = []
        for i in range(length):
            t = traces[i, :].reshape(length, 1)
            x, y = FFT_H_int_all_traces_one_plot(t)
            lower_kvals.append(x)
            lower_Pxx.append(y)
        return lower_kvals, lower_Pxx

    lower_kvals_x, lower_Pxx_x = process_traces(surface, x_length)
    lower_kvals_y, lower_Pxx_y = process_traces(surface.T, y_length)

    comb_px = np.vstack((lower_Pxx_x, lower_Pxx_y))
    Px = np.mean(comb_px, axis=0)
    kvals = lower_kvals_x[0]

    vals_to_remove = 1
    init_rem = 2

    Px_mean_no_first = Px[init_rem:-vals_to_remove]
    kvals_mean_no_last = kvals[:-(vals_to_remove + init_rem)]
    
    coefficients=np.polyfit(np.log10(kvals_mean_no_last),np.log10(Px_mean_no_first),1) 
    polynomial=np.poly1d(coefficients)
    log10_x_fit=polynomial(np.log10(kvals_mean_no_last))
    
    log_kvals = np.log10(kvals_mean_no_last)
    log_Px = np.log10(Px_mean_no_first)
    res = stats.linregress(log_kvals, log_Px)

    slope = res.slope
    c = 10**res.intercept
    H_out = (slope + 1) / -2
    K = kvals_mean_no_last
    N = (len(kvals_mean_no_last) * 2)
    a = ((2 * np.sqrt(2)) / N) * np.sqrt(c)
    N = (N / 2) - 1
    b = np.sqrt(np.sum((K**-(H_out + 0.5)) * np.sin(np.pi * (K / N)))**2)
    std_h = a * b

    plt.loglog(kvals_mean_no_last,10**log10_x_fit,color='r', label='best fit')
    # plt.plot(x,y)
    plt.plot([],[],' ',label='slope = %s'% float('%.3g' % slope))
    plt.plot([],[],' ',label='H = %s'% float('%.3g' % H_out))
    plt.plot([],[],' ',label='Sp = %s'% float('%.4g' % std_h))
    # plt.plot([],[],' ',label='intercept = %.2e' % c)
    plt.plot([],[],' ',label='intercept = %s'% float('%.3g' % c))
    plt.plot(kvals_mean_no_last,Px_mean_no_first, color='C0')
    plt.xlabel('wavenumber, k', fontsize=13)
    plt.ylabel('P(k)', fontsize=13)
    # plt.text(2, 4e3, fig_lab,fontsize=14) #a
    plt.text(2, 2e3, fig_lab,fontsize=14) #b
    plt.legend(fontsize=12)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.show()
    
plot = False
if plot == True:
    # plot_unfiltered_PSD(upper,'(a)','calc_H_PSD_unfiltered/upper_mean_psd_no_filter.png')
    plot_unfiltered_PSD(lower,'(b)','calc_H_PSD_unfiltered/lower_mean_psd_no_filter.png')
