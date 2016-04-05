
load series_data_1e-6

for n=50:50:500
    
    Series_atrium = fourier_series_function_handle(a_0_atrium, a_n_atrium, b_n_atrium, n, true_cycle_length); 
    
    Series_ventricle = fourier_series_function_handle(a_0_ventricle, a_n_ventricle, b_n_ventricle, n, true_cycle_length); 
    
    t = 0:.0001:true_cycle_length*2; 
    vals_atrium_series    = Series_atrium(t); 
    vals_ventricle_series = Series_ventricle(t); 
    figure; 
    plot(t, vals_atrium_series, '--'); 
    hold on; 
    plot(t, vals_ventricle_series); 
    legend('atrial pressure', 'ventricular pressure', 'location', 'NorthWest'); 
    title(sprintf('series data n = %d', n)); 
    
end 





