

% Taken from beat 1, p. 227, fig 2 
% 'dynamics of left ventricular filling' Edward Yellin 
% In Cardiac Mechanics and Function in the Normal and Diseased Heart 

cycle_length = 0.8; 

points_one_cycle = [0.0,   0; 
0.02, -4; 
0.06, 2; 
0.40, 6; 
0.53, 14; 
0.58, 110; 
0.75, 120; 
cycle_length, 8]; 

dt = 1e-4; 
bump_radius = .05; 
n_fourier_coeffs = 1000; 
plots = false; 

[a_0_ventricle a_n_ventricle b_n_ventricle Series_ventricle] = series_and_smooth(points_one_cycle, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_ventricle(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Ventricular pressure')
xlabel('t')
ylabel('p (mmHg)')
printfig(fig, 'ventricular_pressure_yellin')


% t = 0:dt:(3*cycle_length); 
% vals_ventricle_series = Series_ventricle(t); 
% fig = figure; 
% plot(t, vals_ventricle_series, 'k'); 
% title('Ventricular pressure')
% xlabel('t')
% ylabel('p (mmHg)')
% printfig(fig, 'ventricular_pressure_yellin_three_cycles')
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




points_one_cycle_atrium = [0.0, 22.5; 
0.06, 4; 
0.40, 7; 
0.47, 20; 
0.53, 5; 
0.58, 7; 
0.7,  10; 
cycle_length, 22.5]; 

[a_0_atrium a_n_atrium b_n_atrium Series_atrium] = series_and_smooth(points_one_cycle_atrium, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_atrium_series = Series_atrium(t); 
fig = figure; 
plot(t, vals_atrium_series, 'k'); 
title('Atrial pressure')
xlabel('t')
ylabel('p (mmHg)')
printfig(fig, 'atrial_pressure_yellin')



fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
hold on
plot(t, vals_atrium_series, 'k'); 
title('Atrial pressure')
xlabel('t')
ylabel('p (mmHg)')
printfig(fig, 'both_pressure_yellin')


fig = figure; 
p_diff = vals_atrium_series - vals_ventricle_series; 
plot(t, p_diff , 'k'); 
hold on 
plot(t, 0*p_diff , 'k--'); 
title('Pressure difference')
xlabel('t')
ylabel('p (mmHg)')
printfig(fig, 'pressure_diff_yellin')



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



save('series_data_yellin')











% 
% 
% n_coeffs_to_output = 400; 
% 
% n = n_coeffs_to_output; 
% a_n = a_n(1:n);
% b_n = b_n(1:n);
% series_no_array = @(t) a_0 + sum(a_n .* cos((2*pi/cycle_length) * (1:n) .* t)' + ...  
%                                  b_n .* sin((2*pi/cycle_length) * (1:n) .* t)' );   
% 
% Series_truncated = @(t) arrayfun(series_no_array, t); 
% 
% t = 0:dt:cycle_length; 
% vals_ventricle_series = Series_truncated(t); 
% fig = figure; 
% plot(t, vals_ventricle_series, 'k'); 
% title('Ventricular pressure, truncated series')
% xlabel('t')
% ylabel('p (mmHg)')
% 
% 
% output_series_coeffs_to_txt(a_0, a_n, b_n, n_coeffs_to_output, cycle_length, file_name); 
% 
% 












