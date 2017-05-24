function [a_0 a_n b_n Series] = series_and_smooth(points_one_cycle, dt, bump_radius, n_fourier_coeffs, plots)
% 
% Takes data, computes piecewise linear interpolant, 
% smooths with convolution with cosine squared bump,
% computes and outputs Fourier coefficients and function handle for series 
% 

debug = false; 


cycle_length = points_one_cycle(end, 1); 

times = (0:dt:(cycle_length-dt))'; 
n_times = length(times); 

% piecewise linear interpolant 
vals = interp1(points_one_cycle(:,1), points_one_cycle(:,2), times); 

if plots 
    fig = figure; 
    plot(times, vals); 
    title('Piecewise linear interpolant'); 
end 

times_three_cycle = [times; times+cycle_length; times+(2*cycle_length)]; 
vals_three_cycle = [vals; vals; vals]; 

% fig = figure; 
% plot(times_three_cycle, vals_three_cycle)
% title('three cycle')

% scaling factor that gives appropriate radius 
h = (2/pi) * bump_radius; 

% this should integrate to one 
cos_bump = @(x) (abs(x) <= bump_radius) .* (1/h) .* (2/pi) .* cos(x/h).^2; 

% want this mesh to be aligned with the previous mesh
% shift by scalar to be approx centered at zero 
mesh_bump = times_three_cycle - times_three_cycle(length(times_three_cycle)/2); 
bump_vals = cos_bump(mesh_bump); 

if debug 
    fig = figure; 
    plot(mesh_bump, bump_vals); 
    title('bump'); 

    approx_integral = dt * sum(bump_vals)  
    integral_by_quad = quad(cos_bump, -1,1) 
end 


smoothed = dt*conv(vals_three_cycle, bump_vals, 'same'); 

smoothed_one_cycle = smoothed( (n_times+1) : (2*n_times)); 

if plots
    fig = figure; 
    plot(times, smoothed_one_cycle ); 
    title('After convolution, one cycle')
end 

[a_0 a_n b_n Series] = fourier_series_uniform(times, smoothed_one_cycle, cycle_length, n_fourier_coeffs, dt); 


