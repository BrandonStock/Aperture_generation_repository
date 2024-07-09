import numpy as np
from scipy import signal, stats
import matplotlib.pyplot as plt
import os
import sys
func_path = os.path.join(os.path.dirname(__file__), '../required_functions')
sys.path.append(func_path)
from req_functions import psd_data

##path to pre processed surface data
path_to_surface_data='../Pre_processed_data'
## Load surface data
upper=np.load(path_to_surface_data+'/upper.npy')
lower=np.load(path_to_surface_data+'/lower.npy')

def calc_H_psd_full_surf(surface):
    kvals_lower,abins_lower=psd_data(surface)
    vals_to_remove=1
    x=kvals_lower[:-(vals_to_remove+1)]
    y=abins_lower[1:-vals_to_remove]

    log_kvals = np.log10(x)
    log_Px = np.log10(y)
    res = stats.linregress(log_kvals, log_Px)

    slope = res.slope
    c = 10**res.intercept
    H_out = (slope + 1) / -2
    K = x
    N = (len(x) * 2)**2
    a = ((2 * np.sqrt(2)) / N) * np.sqrt(c)
    N = (N / 2) - 1
    b = np.sqrt(np.sum((K**-(H_out + 0.5)) * np.sin(np.pi * (K / N)))**2)
    std_h = a * b
    return H_out, std_h
    
## lower
H_lower, Sp_lower=calc_H_psd_full_surf(lower)
print (H_lower, Sp_lower)
np.save('calc_H_PSD_full_surf_unfiltered/H_lower.npy',H_lower)
np.save('calc_H_PSD_full_surf_unfiltered/Sp_lower.npy',Sp_lower)

## upper
H_upper, Sp_upper=calc_H_psd_full_surf(upper)
print (H_upper, Sp_upper)
np.save('calc_H_PSD_full_surf_unfiltered/H_upper.npy',H_upper)
np.save('calc_H_PSD_full_surf_unfiltered/Sp_upper.npy',Sp_upper)


def plot_unfiltered_full_surf_PSD(surface, fig_lab, output):
    kvals_lower,abins_lower=psd_data(surface)
    vals_to_remove=1
    x=kvals_lower[:-(vals_to_remove+1)]
    y=abins_lower[1:-vals_to_remove]
    
    coefficients=np.polyfit(np.log10(x),np.log10(y),1)
    polynomial=np.poly1d(coefficients)
    log10_x_fit=polynomial(np.log10(x))    
    slope=coefficients[0]
    c=10**coefficients[1]
    
    H_out1=(slope+1)/-2
    N=(len(x)*2)**2
    K=x[:-1]
    a=((2*np.sqrt(2))/N) * np.sqrt(c)
    N=(N/2)-1 
    b=np.sqrt(np.sum((K**-(H_out1+1/2))*(np.sin(np.pi*(K/N))))**2)
    std_h1=a*b
    
    plt.loglog(x,10**log10_x_fit,'r', label='best fit')
    plt.plot(x,y)
    plt.plot([],[],' ',label='slope = %s'% float('%.3g' % slope))
    plt.plot([],[],' ',label='H = %s'% float('%.3g' % H_out1))
    plt.plot([],[],' ',label='Sp = %s'% float('%.4g' % std_h1))
    plt.plot([],[],' ',label='intercept = %.2e' % c)
    # plt.plot([],[],' ',label='intercept = %s'% float('%.3g' % c))
    plt.xlabel('wavenumber, k', fontsize=13)
    plt.ylabel('P(k)', fontsize=13)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.text(3, 9e12, fig_lab,fontsize=14)
    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.show()
    
plot = False
if plot == True:
    plot_unfiltered_full_surf_PSD(lower,'(b)','calc_H_PSD_full_surf_unfiltered/lower_full_psd_no_filter.png')        
    plot_unfiltered_full_surf_PSD(upper,'(a)','calc_H_PSD_full_surf_unfiltered/upper_full_psd_no_filter.png')        

