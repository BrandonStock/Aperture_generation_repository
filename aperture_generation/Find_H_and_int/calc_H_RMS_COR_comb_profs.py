import numpy as np
import matplotlib.pyplot as plt
# from scipy.stats import iqr
from scipy import signal
import os
import sys
func_path = os.path.join(os.path.dirname(__file__), '../required_functions')
sys.path.append(func_path)
from req_functions import RMS_COR

def calculate_RMS_COR_comb_prof(surface, plot=None):

    x_length = surface.shape[1]
    y_length = surface.shape[0]

    # Prepare data for processing
    Tx = [surface[i, :].reshape(x_length, 1) for i in range(y_length)]
    Ty = [surface[:, i].reshape(y_length, 1) for i in range(x_length)]
    all_t = Tx + Ty

    std_full = []
    for trace in all_t:
        # Detrend the data
        trace = signal.detrend(trace, axis=0)

        nofData = len(trace)
        maxBin = int(nofData / 2)

        # Generate step sizes
        steps = np.arange(1, int(0.1 * nofData))  # Using 10% of max length scale

        # Calculate the std of height differences for different step sizes
        std = []
        for step in steps:
            zero_array = np.zeros((step, 1))
            offset_value = np.vstack([zero_array, trace])[:-step]
            diff = trace[step:] - offset_value[step:]
            std.append(np.std(diff))

        std_full.append(std)

    # Flatten std_full list
    std = np.hstack(std_full)

    # Calculate slope and intercept
    new_steps = np.hstack([steps] * len(all_t))
    m, c = np.polyfit(np.log(new_steps), np.log(std), 1)

    # Correction based on Marsch 2021 paper
    if plot == True:
        return new_steps,std,m,c
    
    if m > 0.5:
        m = np.log(m) + 1.18
    c = np.exp(c)

    return m, c


def plot_results(new_steps, std, m, c, letter, output_file):
    """
    Plot the results of the scaling parameter calculation.

    Args:
        new_steps (numpy.ndarray): Array of step sizes.
        std (numpy.ndarray): Array of standard deviations.
        m (float): Slope.
        c (float): Intercept.
        output_file (str): Path to save the plot.
    """
    y_fit = np.exp(m * np.log(new_steps) + c)
    if m > 0.5:
        m = np.log(m) + 1.18
    c = np.exp(c)
       
    plt.scatter(new_steps, std, label='Data', color='red', marker='+')
    plt.plot(new_steps, y_fit, label='Inferred slope')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Length difference',fontsize=12)
    plt.ylabel('Std of height difference',fontsize=12)
    plt.plot([], [], ' ', label='Slope = {:.3g}'.format(m))
    plt.plot([], [], ' ', label='Intercept = {:.3g}'.format(c))
    plt.text(15, 2, letter,fontsize=14)
    plt.legend(fontsize=12)
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.show()
    

    
##path to pre processed surface data
path_to_surface_data='../Pre_processed_data'
## Load surface data
upper=np.load(path_to_surface_data+'/upper.npy')
lower=np.load(path_to_surface_data+'/lower.npy')

## calculate H and Sp lower
H_lower,Sp_lower=calculate_RMS_COR_comb_prof(lower)
print (H_lower,Sp_lower)
np.save('calc_H_RMS_COR_comb_profs/H_lower.npy',H_lower)
np.save('calc_H_RMS_COR_comb_profs/Sp_lower.npy',Sp_lower)

## calculate H and Sp upper
H_upper,Sp_upper=calculate_RMS_COR_comb_prof(upper)
print (H_upper,Sp_upper)
np.save('calc_H_RMS_COR_comb_profs/H_upper.npy',H_upper)
np.save('calc_H_RMS_COR_comb_profs/Sp_upper.npy',Sp_upper)


make_plots=False
if make_plots==True:
    new_steps,std,m,c=calculate_RMS_COR_comb_prof(lower, plot=True)
    lower_plot=plot_results(new_steps,std,m,c,'(b)','calc_H_RMS_COR_comb_profs/rms_cor_lower_allXY_profile.png')

    new_steps,std,m,c=calculate_RMS_COR_comb_prof(upper, plot=True)
    upper_plot=plot_results(new_steps,std,m,c,'(a)','calc_H_RMS_COR_comb_profs/rms_cor_upper_allXY_profile.png')


