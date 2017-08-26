function ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix)



file_name = strcat(base_name, suffix, '.txt'); 

bump_radius = .05; 
n_fourier_coeffs = 1000; 
plots = false; 

[a_0_ventricle a_n_ventricle b_n_ventricle Series_ventricle] = series_and_smooth(points_one_cycle_ventricle, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_ventricle(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Ventricular pressure')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('ventricular_pressure_yellin', suffix))

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



[a_0_atrium a_n_atrium b_n_atrium Series_atrium] = series_and_smooth(points_one_cycle_atrium, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_atrium_series = Series_atrium(t); 
fig = figure; 
plot(t, vals_atrium_series, 'k'); 
title('Atrial pressure')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('atrial_pressure_yellin', suffix))



fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
hold on
plot(t, vals_atrium_series, 'k'); 
title('Atrial pressure')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('both_pressure_yellin', suffix))


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
axis([0 2.4 -10 140])
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
legend('Left Ventricle', 'Left Atrium', 'Location', 'NorthWest')
printfig(fig, 'both_pressure_yellin_three_cycles')




fig = figure; 
p_diff = vals_atrium_series - vals_ventricle_series; 
plot(t, p_diff , 'k'); 
hold on 
plot(t, 0*p_diff , 'k--'); 
title('Pressure difference')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('pressure_diff_yellin', suffix))



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



save(strcat('series_data_yellin', suffix)); 


n_coeffs_to_output = 600; 

n = n_coeffs_to_output; 

a_0 = a_0_atrium - a_0_ventricle; 
a_n = a_n_atrium - a_n_ventricle; 
b_n = b_n_atrium - b_n_ventricle; 

a_n = a_n(1:n);
b_n = b_n(1:n);
series_no_array = @(t) a_0 + sum(a_n .* cos((2*pi/cycle_length) * (1:n) .* t)' + ...  
                                 b_n .* sin((2*pi/cycle_length) * (1:n) .* t)' );   

Series_truncated = @(t) arrayfun(series_no_array, t); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_truncated(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Pressure diff, truncated series')
xlabel('t')
ylabel('p (mmHg)')

fig = figure; 
semilogy( abs(a_n), 'k')
hold on 
semilogy( abs(b_n), ':k')
legend('atrium cos', 'atrium sin')
xlabel('n')
ylabel('|a_n|, |b_n|')
title('Modulus of Fourier coefficients')


output_series_coeffs_to_txt(a_0, a_n, b_n, n_coeffs_to_output, cycle_length, file_name); 














