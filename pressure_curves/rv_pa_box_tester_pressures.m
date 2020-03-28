
% quadrature spacing 
debug = false; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 
plots = true; 

table = readtable('HealthyNativePressures.xlsx');

times = table.Time; 
rv_pressure = table.RightVentriclePressure_Inlet_; 
pa_pressure = table.MainPulmonaryArteryPressure_Outlets_; 

rv_pressure = rv_pressure + 10; 

% periodic wrap 
times = [0; times]; 
rv_pressure = [rv_pressure(end); rv_pressure]; 
pa_pressure = [pa_pressure(end); pa_pressure]; 

cycle_length = max(times); 
base_name = 'fourier_coeffs';
bump_radius = .017; 
n_fourier_coeffs = 600; 

points_one_cycle_right_ventricle = [times, rv_pressure]; 
points_one_cycle_pa              = [times, pa_pressure]; 

manual_style = false; 
if manual_style 
    cycle_length = 0.8; 
    base_name = 'fourier_coeffs';
    bump_radius = .05; 
    n_fourier_coeffs = 600; 

    points_one_cycle_right_ventricle = [0.0,   0; 
    0.02, -4; 
    0.06, 2; 
    0.40, 3; % 6 ; 
    0.53, 14  / 3.8235; 
    0.58, 120 / 3.8235; 
    0.75, 130 / 3.8235; 
    cycle_length, 8]; 

    points_one_cycle_aorta = [0.0,  120;  
    0.569, 80; 
    0.58, 80; 
    0.59, 98; 
    0.60, 110;
    0.61, 117;
    0.72, 122;
    0.73, 122;
    % 0.745, 119; % remove notch 
    % 0.755, 100; 
    0.76, 122;
    cycle_length, 120]; 

    min_p_pa = 8; 
    max_p_pa = 24; 
    pulse_pressure_pa = max_p_pa - min_p_pa; 
    p = points_one_cycle_aorta(:,2); 
    pulse_pressure_aortic = max(p)-min(p); 
    pressure_pa = (p - min(p))*(pulse_pressure_pa/pulse_pressure_aortic) + min_p_pa; 

    points_one_cycle_pa = [points_one_cycle_aorta(:,1) pressure_pa]; 
end 




suffix_right_ventricle = "_right_ventricle"; 

file_name = strcat(base_name, suffix_right_ventricle, '.txt'); 
[a_0_right_ventricle a_n_right_ventricle b_n_right_ventricle Series_right_ventricle] = series_and_smooth(points_one_cycle_right_ventricle, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_right_ventricle_series = Series_right_ventricle(t); 
fig = figure; 
hold on 
plot(t, vals_right_ventricle_series, 'k'); 

output_series_coeffs_to_txt(a_0_right_ventricle, a_n_right_ventricle, b_n_right_ventricle, n_fourier_coeffs, cycle_length, file_name); 

suffix_pa = "_pa"; 

file_name = strcat(base_name, suffix_pa, '.txt'); 
[a_0_pa a_n_pa b_n_pa Series_pa] = series_and_smooth(points_one_cycle_pa, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_pa_series = Series_pa(t); 
figure(fig);
plot(t, vals_pa_series, 'k'); 
title('RV PA pressure')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, 'rv_pa_pressure')

output_series_coeffs_to_txt(a_0_pa, a_n_pa, b_n_pa, n_fourier_coeffs, cycle_length, file_name); 

fig = figure; 
plot(t, vals_right_ventricle_series - vals_pa_series, 'k')
title('Pressure differences (mmHg)')



