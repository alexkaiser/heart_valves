% Taken from beat 1, p. 227, fig 2 
% 'dynamics of left ventricular filling' Edward Yellin 
% In Cardiac Mechanics and Function in the Normal and Diseased Heart 




cycle_length = 0.8; 

% quadrature spacing 
debug = false; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 

plots = false; 

base_name = 'fourier_coeffs';


points_one_cycle_ventricle = [0.0,   0; 
0.02, -4; 
0.06, 2; 
0.40, 6; 
0.53, 14; 
0.58, 120 + 10; 
0.75, 130 + 10; 
cycle_length, 8]; 

points_one_cycle_atrium = [0.0, 24.555; 
0.06, 4; 
0.40, 7; 
0.47, 20; 
0.53, 5; 
0.58, 7; 
0.7,  10; 
cycle_length, 24.555]; 

points_one_cycle_aorta = [0.0,  120;  
0.569, 80; 
0.58, 80; 
0.59, 98; 
0.60, 106;
0.61, 115;
0.72, 124;
0.73, 122;
0.745, 119; 
0.755, 100; 
0.76, 122;
cycle_length, 120]; 

suffix = ''; 


file_name = strcat(base_name, suffix, '.txt'); 

bump_radius = .05; 
n_fourier_coeffs = 600; 
% plots = false; 

[a_0_ventricle a_n_ventricle b_n_ventricle Series_ventricle] = series_and_smooth(points_one_cycle_ventricle, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_ventricle(t); 
% fig = figure; 
% plot(t, vals_ventricle_series, 'k'); 
% title('Ventricular pressure')
% xlabel('t')
% ylabel('p (mmHg)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, strcat('ventricular_pressure_yellin', suffix))

% 
% fig = figure; 
% semilogy( abs(a_n_ventricle), 'k')
% hold on 
% semilogy( abs(b_n_ventricle), ':k')
% legend('ventricle cos', 'ventricle sin')
% xlabel('n')
% ylabel('|a_n|, |b_n|')
% % title('Modulus of Fourier coefficients')
% printfig(fig, 'coefficients')


% n_points_temp = 101; 
% times_temp = linspace(0,cycle_length,n_points_temp); 
% vals_ventricle = Series_ventricle(times_temp); 
% f = fopen('pressure_vals_coarse.txt', 'w'); 
% for i=1:n_points_temp 
%     fprintf(f, '%f, %f\n', times_temp(i), vals_ventricle(i)); 
% end 
% 
% n_points_temp = 10001; 
% times_temp = linspace(0,cycle_length,n_points_temp); 
% vals_ventricle = Series_ventricle(times_temp); 
% f = fopen('pressure_vals_fine.txt', 'w'); 
% for i=1:n_points_temp 
%     fprintf(f, '%f, %f\n', times_temp(i), vals_ventricle(i)); 
% end 




[a_0_atrium a_n_atrium b_n_atrium Series_atrium] = series_and_smooth(points_one_cycle_atrium, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_atrium_series = Series_atrium(t); 
% fig = figure; 
% plot(t, vals_atrium_series, 'k'); 
% title('Atrial pressure')
% xlabel('t')
% ylabel('p (mmHg)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, strcat('atrial_pressure_yellin', suffix))

bump_radius = .015; 

[a_0_aorta a_n_aorta b_n_aorta Series_aorta] = series_and_smooth(points_one_cycle_aorta, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_aorta_series = Series_aorta(t); 
% fig = figure; 
% plot(t, vals_aorta_series, 'k'); 
% title('Aorta pressure')
% xlabel('t')
% ylabel('p (mmHg)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, strcat('aorta_pressure_yellin', suffix))

min_pressure_aorta = min(vals_aorta_series)
max_pressure_aorta = max(vals_aorta_series)
mean_pressure_aorta = mean(vals_aorta_series)

time_zero_pressure_aorta = vals_aorta_series(1)

vals_plus_one  = [vals_aorta_series(2:end), vals_aorta_series(1)]; 
vals_minus_one = [vals_aorta_series(end), vals_aorta_series(1:(end-1))]; 
dp_dt = (vals_plus_one - vals_minus_one)/(2*dt); 

min_dp_dt_aorta = min(dp_dt)
max_dp_dt_aorta = max(dp_dt)

fig = figure; 
plot(t, dp_dt, 'k'); 
hold on
title('dP/dt')
xlabel('t')
ylabel('dp/dt (mmHg/s)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('dp_dt', suffix))

fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
hold on
plot(t, vals_atrium_series, '--k'); 
plot(t, vals_aorta_series, ':k'); 
title('Three pressures')
xlabel('t')
ylabel('p (mmHg)')
legend('LV', 'LA', 'Aorta','Location', 'SouthEast')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('both_pressure_yellin', suffix))


fig = figure; 
plot(t, vals_ventricle_series - vals_aorta_series, 'k'); 
hold on
plot(t, 0*vals_ventricle_series, 'k'); 
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
title('LV minus aorta')
printfig(fig, 'pressure_difference_lv_aorta')


t = 0:dt:(3*cycle_length); 
vals_ventricle_series = Series_ventricle(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Pressures')
xlabel('t')
ylabel('p (mmHg)')
hold on 
vals_atrium_series = Series_atrium(t); 
plot(t, vals_atrium_series, '--k');
vals_aorta_series = Series_aorta(t); 
plot(t, vals_aorta_series, ':k'); 
axis([0 2.4 -10 140])
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
legend('Left Ventricle', 'Left Atrium', 'Location', 'NorthWest')
printfig(fig, 'both_pressure_yellin_three_cycles')



% fig = figure; 
% p_diff = vals_atrium_series - vals_ventricle_series; 
% plot(t, p_diff , 'k'); 
% hold on 
% plot(t, 0*p_diff , 'k--'); 
% title('Pressure difference')
% xlabel('t')
% ylabel('p (mmHg)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% printfig(fig, strcat('pressure_diff_yellin', suffix))



% t = 0:dt:(3*cycle_length); 
% vals_atrium_series = Series_atrium(t); 
% fig = figure; 
% plot(t, vals_atrium_series, 'k'); 
% title('Atrial pressure')
% xlabel('t')
% ylabel('p (mmHg)')
% printfig(fig, 'atrial_pressure_yellin_three_cycles')
% 
% fig = figure; 
% semilogy( abs(a_n), 'k')
% hold on 
% semilogy( abs(b_n), ':k')
% legend('atrium cos', 'atrium sin')
% xlabel('n')
% ylabel('|a_n|, |b_n|')
% title('Modulus of Fourier coefficients')


% 
% n_coeffs_to_output = 600; 
% 
% n = n_coeffs_to_output; 
% 
% a_0 = a_0_atrium - a_0_ventricle; 
% a_n = a_n_atrium - a_n_ventricle; 
% b_n = b_n_atrium - b_n_ventricle; 
% 
% a_n = a_n(1:n);
% b_n = b_n(1:n);
% series_no_array = @(t) a_0 + sum(a_n .* cos((2*pi/cycle_length) * (1:n) .* t)' + ...  
%                                  b_n .* sin((2*pi/cycle_length) * (1:n) .* t)' );   
% 
% Series_truncated = @(t) arrayfun(series_no_array, t); 

% t = 0:dt:cycle_length; 
% vals_ventricle_series = Series_truncated(t); 
% fig = figure; 
% plot(t, vals_ventricle_series, 'k'); 
% title('Pressure diff, truncated series')
% xlabel('t')
% ylabel('p (mmHg)')

% fig = figure; 
% semilogy( abs(a_n), 'k')
% hold on 
% semilogy( abs(b_n), ':k')
% legend('atrium cos', 'atrium sin')
% xlabel('n')
% ylabel('|a_n|, |b_n|')
% title('Modulus of Fourier coefficients')


% output_series_coeffs_to_txt(a_0, a_n, b_n, n_coeffs_to_output, cycle_length, file_name); 


n_coeffs_to_output = n_fourier_coeffs; 

output_series_coeffs_to_txt(a_0_ventricle, a_n_ventricle, b_n_ventricle, n_coeffs_to_output, cycle_length, 'fourier_coeffs_ventricle.txt'); 
output_series_coeffs_to_txt(a_0_atrium,    a_n_atrium,    b_n_atrium,    n_coeffs_to_output, cycle_length, 'fourier_coeffs_atrium.txt'); 
output_series_coeffs_to_txt(a_0_aorta,     a_n_aorta,     b_n_aorta,     n_coeffs_to_output, cycle_length, 'fourier_coeffs_aorta.txt'); 


a_0 = a_0_aorta - a_0_ventricle; 
a_n = a_n_aorta - a_n_ventricle; 
b_n = b_n_aorta - b_n_ventricle; 

output_series_coeffs_to_txt(a_0, a_n, b_n, n_coeffs_to_output, cycle_length, 'fourier_coeffs_aorta_minus_ventricle.txt'); 



% close all 
% save(strcat('series_data', suffix)); 






