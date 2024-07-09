## file to define variables

'number of realisations that are generated, total number = number_of_realisations - 1'
number_of_realisations = 101

'dimensions of generation, size defined by 2**N, total width 2*2**N+1'
dimension_N = 8

'number of processors used to execute the number swapping algorithm'
number_of_processors = 9

'order of the polynomial regression over the correlation (1-PSDR), might need changing for other fractures'
poly_deg = 3

'define which method is used to calculate H and Sp'
calc_H_PSD_filtered = False
calc_H_PSD_unfiltered = False
calc_H_PSD_full_surf_filtered = False
calc_H_PSD_full_surf_unfiltered = False
calc_H_RMS_COR = True
calc_H_RMS_COR_comb_profs = False
