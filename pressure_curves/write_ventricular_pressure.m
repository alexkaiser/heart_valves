


load series_data_1e-6.mat

file_name = 'fourier_coeffs_ventricle.txt'; 

output_series_coeffs_to_txt(a_0_ventricle, a_n_ventricle, b_n_ventricle, n, true_cycle_length, file_name); 


t = 0:.0001:true_cycle_length; 
vals_ventricle_series = Series_ventricle(t); 
figure; 
hold on; 
plot(t, vals_ventricle_series); 
title('Ventricular pressure from fourier seriess, one cycle')

t = 0:.0001:true_cycle_length*2; 
vals_ventricle_series = Series_ventricle(t); 
figure; 
hold on; 
plot(t, vals_ventricle_series); 
title('Ventricular pressure from fourier seriess, two cycles')




