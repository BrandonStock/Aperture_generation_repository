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
    
    remove_vals=np.arange(1,900)
    score_lower=[]
    Hs=[]
    Sps=[]
    for i in range (len(remove_vals)):
        x=kvals_lower[:-(remove_vals[i]+2)]
        y=abins_lower[2:-remove_vals[i]]
        res = stats.linregress(np.log10(x),np.log10(y))
        slope=res.slope
        c=10**res.intercept
        H_out=(slope+1)/-2
        K=x 
        N=(len(x)*2)**2
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
H_lower, Sp_lower=calc_H_psd_full_surf(lower)
print (H_lower, Sp_lower)
np.save('calc_H_PSD_full_surf_filtered/H_lower.npy',H_lower)
np.save('calc_H_PSD_full_surf_filtered/Sp_lower.npy',Sp_lower)

## upper
H_upper, Sp_upper=calc_H_psd_full_surf(upper)
print (H_upper, Sp_upper)
np.save('calc_H_PSD_full_surf_filtered/H_upper.npy',H_upper)
np.save('calc_H_PSD_full_surf_filtered/Sp_upper.npy',Sp_upper)

    
def calc_H_psd_full_surf_filtering(surface):
    kvals_lower,abins_lower=psd_data(surface)
    vals_to_remove=1
    x=kvals_lower[:-(vals_to_remove+2)]
    y=abins_lower[2:-vals_to_remove]

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
    
    remove_vals=np.arange(1,900)
    score_lower=[]
    Hs=[]
    Sps=[]
    for i in range (len(remove_vals)):
        x=kvals_lower[:-(remove_vals[i]+2)]
        y=abins_lower[2:-remove_vals[i]]
        res = stats.linregress(np.log10(x),np.log10(y))
        slope=res.slope
        c=10**res.intercept
        H_out=(slope+1)/-2
        K=x 
        N=(len(x)*2)**2
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
        
def plot_filtered_full_surf_PSD(surface,fig_lab, vals_to_remove,output):
    kvals_lower,abins_lower=psd_data(surface)
    # vals_to_remove=1
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
    plt.text(2, 15e12, fig_lab,fontsize=14)
    plt.savefig(output, dpi=600, bbox_inches='tight')
    plt.show()
    

plot_figs = False
if plot_figs == True:
        
    # ## Hs and Sps lower
    Hs_lower, Sps_lower=calc_H_psd_full_surf_filtering(lower)
    
    # ## Hs and Sps upper
    Hs_upper, Sps_upper=calc_H_psd_full_surf_filtering(upper)
    
    ## plot H vs Sp
    plt.scatter(Hs_lower,Sps_lower,label='lower', marker='.')
    plt.scatter(Hs_upper,Sps_upper,label='upper', marker='.')
    plt.ylabel('Sp', fontsize=13)
    plt.xlabel('H',fontsize=13)
    plt.legend(fontsize=12)
    plt.text(0.455, 7.5, '(a)',fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    # plt.yscale('log')
    plt.savefig('calc_H_PSD_full_surf_filtered/Sp_vs_H_full_surf_filtered.png', dpi=600, bbox_inches='tight')
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
    plt.text(0.2, 0.54, '(b)',fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.savefig('calc_H_PSD_full_surf_filtered/H_vs_wavelength_full_surf_filtered.png', dpi=600, bbox_inches='tight')
    plt.show()
    
    
    plot_filtered_full_surf_PSD(lower,'(d)',idx_lower[0][0],'calc_H_PSD_full_surf_filtered/lower_full_psd_filter.png')
    plot_filtered_full_surf_PSD(upper,'(c)',idx_upper[0][0],'calc_H_PSD_full_surf_filtered/upper_full_psd_filter.png')
    
    
    
