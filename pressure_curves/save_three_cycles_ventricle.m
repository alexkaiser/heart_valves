
load series_data_yellin.mat

dt = 1e-3; 

t = 0:dt:(cycle_length*3); 

vals_ventricle_series = Series_ventricle(t); 

save('ventricle_three_cycles', 't', 'vals_ventricle_series'); 


