

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

dt = 1e-5; 
times = (0:dt:(cycle_length-dt))'; 

n_times = length(times); 

vals = interp1(  points_one_cycle(:,1), points_one_cycle(:,2), times); 

fig = figure; 
plot(times, vals); 
title('ventricular pressure, piecewise linear')


times_three_cycle = [times; times+cycle_length; times+(2*cycle_length)]; 

vals_three_cycle = [vals; vals; vals]; 

% fig = figure; 
% plot(times_three_cycle, vals_three_cycle)
% title('three cycle')

% total width of bump 
width = .15; 

% this should integrate to one 
cos_bump = @(x) (abs(x) <= width/4) .* (pi/width) .* cos( (2*pi/width) * x); 

% want this mesh to be aligned with the previous mesh
% shift by scalar to be approx centered at zero 
mesh_bump = times_three_cycle - times_three_cycle(length(times_three_cycle)/2); 
bump_vals = cos_bump(mesh_bump); 

% fig = figure; 
% plot(mesh_bump, bump_vals); 
% title('bump'); 

approx_integral = dt * sum(bump_vals)  

integral_by_quad = quad(cos_bump, -1,1); 


smoothed = dt*conv(vals_three_cycle, bump_vals, 'same'); 

% fig = figure; 
% plot(smoothed); 
% title('after convolution')

smoothed_one_cycle = smoothed( (n_times+1) : (2*n_times)); 

fig = figure; 
plot(times, smoothed_one_cycle ); 
title('after convolution, one cycle')


n_fourier_coeffs = 1000; 
[a_0 a_n b_n Series_ventricle] = fourier_series_uniform(times, smoothed_one_cycle, cycle_length, n_fourier_coeffs, dt); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_ventricle(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Ventricular pressure')
xlabel('t')
ylabel('p (mmHg)')
printfig(fig, 'ventricular_pressure_yellin')


t = 0:dt:(3*cycle_length); 
vals_ventricle_series = Series_ventricle(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Ventricular pressure')
xlabel('t')
ylabel('p (mmHg)')
printfig(fig, 'ventricular_pressure_yellin_three_cycles')

fig = figure; 
semilogy( abs(a_n), 'k')
hold on 
semilogy( abs(b_n), ':k')
legend('ventricle cos', 'ventricle sin')
xlabel('n')
ylabel('|a_n|, |b_n|')
% title('Modulus of Fourier coefficients')
printfig(fig, 'coefficients')

file_name = 'fourier_coeffs_ventricle.txt'; 

save('series_data_yellin')


n_coeffs_to_output = 400; 

n = n_coeffs_to_output; 
a_n = a_n(1:n);
b_n = b_n(1:n);
series_no_array = @(t) a_0 + sum(a_n .* cos((2*pi/cycle_length) * (1:n) .* t)' + ...  
                                 b_n .* sin((2*pi/cycle_length) * (1:n) .* t)' );   

Series_truncated = @(t) arrayfun(series_no_array, t); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_truncated(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Ventricular pressure, truncated series')
xlabel('t')
ylabel('p (mmHg)')


output_series_coeffs_to_txt(a_0, a_n, b_n, n_coeffs_to_output, cycle_length, file_name); 














