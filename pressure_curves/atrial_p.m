
curve_spec_atrium = [ -8.766 -0.108 -16.498 -0.198 -22.344 -0.198 -24.104 0 -31.518 -9.888 -35.226 -9.888 -5.56 0 -66.744 15.452 -90.228 15.452 -27.368 0 -32.136 -13.596 -39.552 -13.596 -8.034 0 -12.362 14.028 -20.238 14.028 -10.044 0 -16.224 -20.826 -37.856 -20.826 -15.138 0 -25.122 11.854 -44.152 14.004 -19.03 2.15 -85.606 0.828 -110.034 0.828 -24.102 0 -31.518 -9.888 -35.226 -9.888 -5.562 0 -66.746 15.452 -90.23 15.452 -27.368 0 -32.136 -13.596 -39.552 -13.596 -8.034 0 -12.36 14.028 -20.238 14.028 -10.044 0 -16.224 -20.826 -37.856 -20.826 -18.31 0 -40.946 18.362 -62.576 12.798];


x_0 = 842.33797; 
y_0 = 419.006; 

% two beats in the current diagram 
min_x =  0.0; 
max_x =  1.6; 
min_y =  1.0; 
max_y =  12.0; 

dt = 1.0e-6; 
% each unit length (in paramter space) is evaluated at this many points 
N_per_unit_length = ceil(1/dt); 


% [x_pressure_atrium y_pressure_atrium] = svg_curve_to_points(curve_spec_atrium, x_0, y_0, min_x, max_x, min_y, max_y); 
[x_pressure_atrium y_pressure_atrium] = eval_bezier_curve_on_path(N_per_unit_length, curve_spec_atrium, x_0, y_0, min_x, max_x, min_y, max_y); 

% figure;
% plot(x_pressure_atrium, y_pressure_atrium)
% title('atrium, no smoothing')

% x_spline = 0:.001:1.6; 
% y_spline_atrium = spline(x_pressure_atrium, y_pressure_atrium, x_spline); 
% 
% figure; 
% plot(x_spline, y_spline_atrium); 
% title('atrial pressure, with spline')


min_x =  0.0; 
max_x =  1.6; 
min_y =  0.0; 
max_y =  65.0;  % reduced here, 20 too low, split difference  
%max_y = 110.0; % physiological  


% probably not the right curve, 
% this is wrong accordinag to parser 
% x_0 = 842.33797;  
% y_0 = 262.864; 
% curve_spec_ventricle = [-26.616 -6.758 -47.618 -13.22 -54.48 -17.098 -11.31 -6.394 -11.124 12.36 -19.776 11.742 -8.652 -0.618 -32.906 -39.732 -56.238 -37.08 -27.194 3.09 -44.498 66.126 -56.858 70.454 -9.95 3.482 -58.39 -1.65 -96.092 -8.284 -37.702 -6.634 -134.11 -28.798 -148.324 -36.832 -11.312 -6.394 -11.124 12.36 -19.776 11.742 -8.652 -0.618 -32.904 -39.732 -56.238 -37.08 -27.192 3.09 -44.498 66.126 -56.856 70.454 -12.358 4.328 -84.826 -3.772 -120.67 -14.278];

x_0 = 842.33797;
y_0 = 432.988; 
curve_spec_ventricle = [-32.566 -0.228 -48.238 -4.746 -52.626 -12.326 -13.488 -23.3 -7.064 -117.176 -15.472 -149.298 -4.922 -18.8 -34.346 -66.154 -70.43 -55.88 -64.85 18.462 -62.418 207.032 -69.06 209.32 -8.884 3.06 -14.372 -13.95 -28.584 -13.95 -25.34 0 -38.86 15.898 -71.528 16.738 -32.668 0.84 -115.414 4.812 -122.212 -6.93 -13.49 -23.3 -7.066 -117.176 -15.474 -149.298 -4.92 -18.8 -34.346 -66.154 -70.43 -55.88 -64.85 18.462 -62.418 207.032 -69.06 209.32 -8.886 3.06 -14.372 -13.95 -28.584 -13.95 -25.338 0 -40.17 16.738 -71.846 16.738]; 

[x_pressure_ventricle y_pressure_ventricle] = eval_bezier_curve_on_path(N_per_unit_length, curve_spec_ventricle, x_0, y_0, min_x, max_x, min_y, max_y); 

% figure; 
% plot(x_pressure_ventricle, y_pressure_ventricle); 
% title('ventricle')

% y_spline_ventricle = spline(x_pressure_ventricle, y_pressure_ventricle, x_spline); 
% 
% figure; 
% plot(x_spline, y_spline_ventricle); 
% title('ventricle pressure, with spline')

figure; 
plot(x_pressure_atrium, y_pressure_atrium, '--'); 
hold on 
plot(x_pressure_ventricle, y_pressure_ventricle); 
legend('atrial pressure', 'ventricular pressure'); 


% manually selected interval, starts at zero p in early diastole  
min_t = 0.5873; 
max_t = 1.4682; 
true_cycle_length = 0.8; 
current_cycle_length = max_t - min_t; 

% traverse in reverse order since the curves are right to left 
count = 1;
x_pressure_atrium_one_cycle = zeros(size(x_pressure_atrium)); 
y_pressure_atrium_one_cycle = zeros(size(y_pressure_atrium)); 

for i = length(x_pressure_atrium):-1:1
    
    if (min_t <= x_pressure_atrium(i)) && (x_pressure_atrium(i) < max_t)
        x_pressure_atrium_one_cycle(count) = x_pressure_atrium(i); 
        y_pressure_atrium_one_cycle(count) = y_pressure_atrium(i); 
        
        count = count + 1; 
    end 
    
end 

x_pressure_atrium_one_cycle = x_pressure_atrium_one_cycle(1:count-1); 
y_pressure_atrium_one_cycle = y_pressure_atrium_one_cycle(1:count-1); 

% scale to .8 s 
x_pressure_atrium_one_cycle = x_pressure_atrium_one_cycle - min(x_pressure_atrium_one_cycle);
x_pressure_atrium_one_cycle = (true_cycle_length / current_cycle_length) * x_pressure_atrium_one_cycle;  


% traverse in reverse order since the curves are right to left 
count = 1;
x_pressure_ventricle_one_cycle = zeros(size(x_pressure_ventricle)); 
y_pressure_ventricle_one_cycle = zeros(size(y_pressure_ventricle)); 

for i = length(x_pressure_ventricle):-1:1
    
    if (min_t <= x_pressure_ventricle(i)) && (x_pressure_ventricle(i) < max_t)
        x_pressure_ventricle_one_cycle(count) = x_pressure_ventricle(i); 
        y_pressure_ventricle_one_cycle(count) = y_pressure_ventricle(i); 
        
        count = count + 1; 
    end 
    
end 

x_pressure_ventricle_one_cycle = x_pressure_ventricle_one_cycle(1:count-1); 
y_pressure_ventricle_one_cycle = y_pressure_ventricle_one_cycle(1:count-1); 

% scale to .8 s 
x_pressure_ventricle_one_cycle = x_pressure_ventricle_one_cycle - min(x_pressure_ventricle_one_cycle);
x_pressure_ventricle_one_cycle = (true_cycle_length / current_cycle_length) * x_pressure_ventricle_one_cycle;  



figure; 
plot(x_pressure_atrium_one_cycle, y_pressure_atrium_one_cycle , '--'); 
hold on 
plot(x_pressure_ventricle_one_cycle , y_pressure_ventricle_one_cycle); 
legend('atrial pressure', 'ventricular pressure', 'location', 'NorthWest'); 
title('single cycle starting at early diastole, zero pressure difference')


% compute the series 
n = 1000; 
[a_0_atrium    a_n_atrium    b_n_atrium    Series_atrium]    = fourier_series(x_pressure_atrium_one_cycle, y_pressure_atrium_one_cycle, true_cycle_length, n); 
[a_0_ventricle a_n_ventricle b_n_ventricle Series_ventricle] = fourier_series(x_pressure_ventricle_one_cycle, y_pressure_ventricle_one_cycle, true_cycle_length, n); 

% files are weirdly large, remove some variables 
clear x_pressure*
clear y_pressure*

% save for fun 
save(sprintf('series_data_ventricular_max_%f_dt_%f.mat', max_y, dt)); 



t = 0:.0001:true_cycle_length; 
vals_atrium_series    = Series_atrium(t); 
vals_ventricle_series = Series_ventricle(t); 
fig = figure; 
plot(t, vals_atrium_series, '--k'); 
hold on; 
plot(t, vals_ventricle_series, 'k'); 
legend('atrial pressure', 'ventricular pressure', 'location', 'NorthWest'); 
title('Driving pressures')
xlabel('t')
ylabel('p (mmHg)')
printfig(fig, 'pressure_curves')

fig = figure; 
plot(t, vals_atrium_series - vals_ventricle_series); 
title('Pressure difference, atrium positive')

t = 0:.0001:true_cycle_length*2; 
vals_atrium_series    = Series_atrium(t); 
vals_ventricle_series = Series_ventricle(t); 
figure; 
plot(t, vals_atrium_series, '--k'); 
hold on; 
plot(t, vals_ventricle_series, 'k'); 
legend('atrial pressure', 'ventricular pressure', 'location', 'NorthWest'); 
title('Driving pressures')
xlabel('t')
ylabel('p (mmHg)')
printfig(fig, 'pressure_curves_two')

fig = figure; 
semilogy( abs(a_n_atrium), '--k' )
hold on 
semilogy( abs(b_n_atrium), '-.k' )
semilogy( abs(a_n_ventricle), 'k')
semilogy( abs(b_n_ventricle), ':k')
legend('atrium cos', 'atrium sin', 'ventricle cos', 'ventricle sin')
xlabel('n')
ylabel('|a_n|, |b_n|')
% title('Modulus of Fourier coefficients')
printfig(fig, 'coefficients')

n_output = 200; 
output_series_to_parser_format(a_0_atrium, a_n_atrium, b_n_atrium, a_0_ventricle, a_n_ventricle, b_n_ventricle, n_output, true_cycle_length); 

a_0 = a_0_atrium - a_0_ventricle; 
a_n = a_n_atrium - a_n_ventricle; 
b_n = b_n_atrium - b_n_ventricle;
file_name = 'fourier_coeffs.txt'; 

output_series_coeffs_to_txt(a_0, a_n, b_n, n, true_cycle_length); 


