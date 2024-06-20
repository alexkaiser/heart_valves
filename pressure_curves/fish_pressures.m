% Taken from beat 1, p. 227, fig 2 
% 'dynamics of left ventricular filling' Edward Yellin 
% In Cardiac Mechanics and Function in the Normal and Diseased Heart 




% cycle_length = 0.8; 

% quadrature spacing 
debug = false; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 

plots = true; 

base_name = 'fourier_coeffs_fish';

% skip header with offsets of 1 
points_one_cycle_ventricle = csvread('Anatomical_Record_ventricular_one_cycle.csv', 1, 0); 

% adjust time to start at zero 
points_one_cycle_ventricle(:,1) = points_one_cycle_ventricle(:,1) - points_one_cycle_ventricle(1,1); 

% cycle length from ventricular pressure 
cycle_length = points_one_cycle_ventricle(end,1) 


stroke_volume_nanol = 266.00; 
NANOL_TO_ML = 1e-6; 

stroke_volume = stroke_volume_nanol * NANOL_TO_ML; 

Q_mean_ml_per_s = stroke_volume / cycle_length

points_one_cycle_aorta = csvread('Anatomical_Record_aortic.csv',1,0); 

points_one_cycle_aorta(:,1) = points_one_cycle_aorta(:,1) - points_one_cycle_aorta(1,1); 

% points_one_cycle_ventricle = [0.0,   0; 
% 0.02, 0; 
% 0.06, 2; 
% 0.40, 6; 
% 0.53, 14; 
% 0.58, 120; 
% 0.75, 130; 
% cycle_length, 8]; 
% 
% points_one_cycle_aorta = [0.0,  120;  
% 0.569, 80; 
% 0.58, 80; 
% 0.59, 98; 
% 0.60, 106;
% 0.61, 115;
% 0.72, 124;
% 0.73, 122;
% 0.745, 119; 
% 0.755, 100; 
% 0.76, 122;
% cycle_length, 120]; 

suffix = ''; 


file_name = strcat(base_name, suffix, '.txt'); 

bump_radius = .05 / 2; 
n_fourier_coeffs = 600; 


[a_0_ventricle a_n_ventricle b_n_ventricle Series_ventricle] = series_and_smooth(points_one_cycle_ventricle, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_ventricle(t); 


bump_radius = .015; 

[a_0_aorta a_n_aorta b_n_aorta Series_aorta] = series_and_smooth(points_one_cycle_aorta, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_aorta_series = Series_aorta(t); 
vals_ventricle_series = Series_ventricle(t); 

min_pressure_aorta = min(vals_aorta_series)
max_pressure_aorta = max(vals_aorta_series)
mean_pressure_aorta = mean(vals_aorta_series)

time_zero_pressure_aorta = vals_aorta_series(1)

vals_plus_one  = [vals_aorta_series(2:end), vals_aorta_series(1)]; 
vals_minus_one = [vals_aorta_series(end), vals_aorta_series(1:(end-1))]; 
dp_dt = (vals_plus_one - vals_minus_one)/(2*dt); 

min_dp_dt_aorta = min(dp_dt)
max_dp_dt_aorta = max(dp_dt)

pressure_diff = vals_aorta_series - vals_ventricle_series; 
figure; 
plot(t, pressure_diff); 

pressure_diff_pos = pressure_diff((pressure_diff > 0)); 

min_fwd_pressure = min(pressure_diff_pos)
max_fwd_pressure = max(pressure_diff_pos)
mean_fwd_pressure = mean(pressure_diff_pos)



% fig = figure; 
% plot(t, dp_dt, 'k'); 
% hold on
% title('dP/dt')
% xlabel('t')
% ylabel('dp/dt (mmHg/s)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, strcat('dp_dt', suffix))
% 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
hold on
plot(t, vals_aorta_series, ':k'); 
title('Three pressures')
xlabel('t')
ylabel('p (mmHg)')
legend('LV', 'LA', 'Aorta','Location', 'SouthEast')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, strcat('both_pressure_yellin', suffix))
% 
% 
% fig = figure; 
% plot(t, vals_ventricle_series - vals_aorta_series, 'k'); 
% hold on
% plot(t, 0*vals_ventricle_series, 'k'); 
% xlabel('t')
% ylabel('p (mmHg)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% title('LV minus aorta')
% printfig(fig, 'pressure_difference_lv_aorta')
% 
% 
% t = 0:dt:(3*cycle_length); 
% vals_ventricle_series = Series_ventricle(t); 
% fig = figure; 
% plot(t, vals_ventricle_series, 'k'); 
% title('Pressures')
% xlabel('t')
% ylabel('p (mmHg)')
% hold on 
% vals_atrium_series = Series_atrium(t); 
% plot(t, vals_atrium_series, '--k');
% vals_aorta_series = Series_aorta(t); 
% plot(t, vals_aorta_series, ':k'); 
% axis([0 2.4 -10 140])
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% legend('Left Ventricle', 'Left Atrium', 'Location', 'NorthWest')
% printfig(fig, 'both_pressure_yellin_three_cycles')
% 


n_coeffs_to_output = n_fourier_coeffs; 

output_series_coeffs_to_txt(a_0_ventricle, a_n_ventricle, b_n_ventricle, n_coeffs_to_output, cycle_length, 'fourier_coeffs_ventricle_fish.txt'); 
% output_series_coeffs_to_txt(a_0_aorta,     a_n_aorta,     b_n_aorta,     n_coeffs_to_output, cycle_length, 'fourier_coeffs_aorta.txt'); 
% 






