import numpy as np
import matplotlib.pyplot as plt
# from scipy.stats import iqr
# import pandas as pd
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
ap=upper-lower

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
    
    remove_vals=np.arange(1,900)
    Hs=[]
    Sps=[]
    for i in range (len(remove_vals)):
        x=kvals[:-(remove_vals[i]+2)]
        y=Px[2:-remove_vals[i]]
        res = stats.linregress(np.log10(x),np.log10(y))
        slope=res.slope
        c=10**res.intercept
        H_out=(slope+1)/-2
        K=x 
        N=(len(x)*2)
        a=((2*np.sqrt(2))/N) * np.sqrt(c)
        N=(N/2)-1 
        b=np.sqrt(np.sum((K**-(H_out+1/2))*(np.sin(np.pi*(K/N))))**2)
        std_h=a*b
        
        Hs.append(H_out)
        Sps.append(std_h)
    
    idx = np.where(Hs==np.max(Hs))
    H_out=Hs[idx[0][0]]
    std_h=Sps[idx[0][0]]
        
    return H_out, std_h

## lower
H_lower, Sp_lower=calculate_scaling_PSD(lower)
print (H_lower, Sp_lower)
np.save('calc_H_PSD_filtered/H_lower.npy',H_lower)
np.save('calc_H_PSD_filtered/Sp_lower.npy',Sp_lower)

## upper
H_upper, Sp_upper=calculate_scaling_PSD(upper)
print (H_upper, Sp_upper)
np.save('calc_H_PSD_filtered/H_upper.npy',H_upper)
np.save('calc_H_PSD_filtered/Sp_upper.npy',Sp_upper)


def change_in_H_and_Sps_filtering(surface):
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
    
    remove_vals=np.arange(1,900)
    Hs=[]
    Sps=[]
    for i in range (len(remove_vals)):
        x=kvals[:-(remove_vals[i]+2)]
        y=Px[2:-remove_vals[i]]
        res = stats.linregress(np.log10(x),np.log10(y))
        slope=res.slope
        c=10**res.intercept
        H_out=(slope+1)/-2
        K=x 
        N=(len(x)*2)
        a=((2*np.sqrt(2))/N) * np.sqrt(c)
        N=(N/2)-1 
        b=np.sqrt(np.sum((K**-(H_out+1/2))*(np.sin(np.pi*(K/N))))**2)
        std_h=a*b
        
        Hs.append(H_out)
        Sps.append(std_h)
    
    idx = np.where(Hs==np.max(Hs))
    H_out=Hs[idx[0][0]]
    std_h=Sps[idx[0][0]]
        
    return Hs, Sps

def plot_filtered_mean_PSD(surface,fig_lab, vals_to_remove,output):
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

    # vals_to_remove = 1
    init_rem = 2

    Px_mean_no_first = Px[init_rem:-vals_to_remove]
    kvals_mean_no_last = kvals[:-(vals_to_remove + init_rem)]
    
    # vals_to_remove=idx[0][0]
    ## remove values from array
    Px_mean_no_first=Px[2:-vals_to_remove]
    kvals_mean_no_last=kvals[:-(vals_to_remove+2)]
    coefficients=np.polyfit(np.log10(kvals_mean_no_last),np.log10(Px_mean_no_first),1) 
    res = stats.linregress(np.log10(kvals_mean_no_last),np.log10(Px_mean_no_first))
    polynomial=np.poly1d(coefficients)
    log10_x_fit=polynomial(np.log10(kvals_mean_no_last))
    slope=res.slope
    c=10**res.intercept
    H_out=(slope+1)/-2
    K=kvals_mean_no_last #kvals_no_last , change to kvals_flat_no_last for means #flat_kvals #x[:-1] 
    N=(len(kvals_mean_no_last)*2)# or maybe len(kvals)*2
    a=((2*np.sqrt(2))/N) * np.sqrt(c)
    N=(N/2)-1 
    b=np.sqrt(np.sum((K**-(H_out+1/2))*(np.sin(np.pi*(K/N))))**2)
    std_h=a*b
    
    plt.loglog(kvals_mean_no_last,10**log10_x_fit,color='r', label='best fit')
    # plt.plot(x,y_fit,label='infered slope')
    plt.plot([],[],' ',label='slope = %s'% float('%.3g' % slope))
    plt.plot([],[],' ',label='H = %s'% float('%.3g' % H_out))
    plt.plot([],[],' ',label='intercept = %s'% float('%.3g' % c))
    plt.plot([],[],' ',label='Sp = %s'% float('%.3g' % std_h))
    plt.plot(kvals_mean_no_last,Px_mean_no_first, color='C0')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('wavenumber, k', fontsize=13)
    plt.ylabel('P(k)', fontsize=13)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    # plt.text(2, 8e3, fig_lab,fontsize=14) #c
    plt.text(2, 3e3, fig_lab,fontsize=14) #d
    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.show()

make_plots = False
if make_plots == True:
        
    # ## Hs and Sps lower
    Hs_lower, Sps_lower=change_in_H_and_Sps_filtering(lower)
    
    # ## Hs and Sps upper
    Hs_upper, Sps_upper=change_in_H_and_Sps_filtering(upper)
    
    ## plot H vs Sp
    plt.scatter(Hs_lower,Sps_lower,label='lower', marker='.')
    plt.scatter(Hs_upper,Sps_upper,label='upper', marker='.')
    plt.ylabel('Sp', fontsize=13)
    plt.xlabel('H',fontsize=13)
    plt.legend(fontsize=12)
    plt.text(0.77, 0.57, '(a)',fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    # plt.yscale('log')
    plt.savefig('calc_H_PSD_filtered/Sp_vs_H_mean_psd_filtered.png', dpi=600, bbox_inches='tight')
    plt.show()
    
    ## plot H vs Sp for different wavelength filtering
    remove_vals=np.arange(1,900)
    wavelengths_removed=1/((len(upper)/2)-remove_vals) * len(upper)/10 
    
    # ## max lower and upper
    idx_lower = np.where(Hs_lower==np.max(Hs_lower))
    idx_upper = np.where(Hs_upper==np.max(Hs_upper))
    
    plt.plot(wavelengths_removed,Hs_lower,label='lower')
    plt.scatter(wavelengths_removed[idx_lower[0][0]],Hs_lower[idx_lower[0][0]], marker='x', label='max H lower')
    plt.plot(wavelengths_removed,Hs_upper,label='upper')
    plt.scatter(wavelengths_removed[idx_upper[0][0]],Hs_upper[idx_upper[0][0]], marker='x', label='max H upper')
    plt.xlabel('Wavelength (mm)',fontsize=13)
    plt.ylabel('H',fontsize=13)
    plt.legend(fontsize=12)
    plt.text(0.5, 0.98, '(b)',fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.savefig('calc_H_PSD_filtered/H_vs_wavelength_mean_psd_filtered.png', dpi=600, bbox_inches='tight')
    plt.show()
    
    plot_filtered_mean_PSD(lower,'(d)',idx_lower[0][0],'calc_H_PSD_filtered/lower_mean_psd_filter.png')
    # plot_filtered_mean_PSD(upper,'(c)',idx_upper[0][0],'calc_H_PSD_filtered/upper_mean_psd_filter.png')



