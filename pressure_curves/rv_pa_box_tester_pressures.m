
% quadrature spacing 
debug = false; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 
plots = true; 

MMHG_TO_CGS = 1333.22368; 

R_proximal_experimental = 4.89   
R_distal_experimental   = 46.15  
C_experimental          = 0.0154 
R_total_experimental    = R_proximal_experimental + R_distal_experimental
tau_experimental        = C_experimental * R_distal_experimental

table = readtable('HealthyNativePressures.xlsx');

times = table.Time; 
rv_pressure = table.RightVentriclePressure_Inlet_; 
pa_pressure = table.MainPulmonaryArteryPressure_Outlets_; 

% rv_pressure = rv_pressure + 20; 

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


flows = [-4.679
3.558
34.035
101.834
180.228
226.853
230.468
192.169
123.601
48.604
4.302
-2.023
-0.795
0.197
0.101
-0.209
-1.054
-2.747
-2.837
-3.134]; 

times_flows = linspace(0, cycle_length, length(flows)); 

flows_spline = interp1(times_flows, flows, t, 'spline'); 

figure; 
plot(times_flows, flows)
hold on 
plot(t, flows_spline) 
title('Flow, experimental trace')
xlabel('t (s)')
ylabel('flow (ml/s)')
legend('pointwise', 'spline')

dt_flows = cycle_length / length(flows)
total_flow = sum(flows) * dt_flows % ml 

Q_mean = total_flow / cycle_length 

Q_mean_L_per_min = Q_mean * 60/1000 

ratio_prox_to_distal_resistors = 0.025; % R_proximal_experimental / R_distal_experimental

decay_time = .4 

P_mean = mean(vals_pa_series) * MMHG_TO_CGS; 

P_max = 42 * MMHG_TO_CGS; % interpolating the decay by eye to the middle of the oscillation 
                          % in the experimental trace of PA pressure 
P_min = 31.25 * MMHG_TO_CGS; 




tol = 1e-12; 

Q_mean_each = Q_mean / 2; 

% total resistance is determined by mean pressure and mean flow 
R_total = P_mean / Q_mean_each;  

% ratio of resistors is constant 
% resistors sum to total resistance 
R_distal = R_total / (1.0 + ratio_prox_to_distal_resistors); 
R_proximal = R_total - R_distal; 

if abs(R_distal + R_proximal - R_total) > tol 
    error('resistors not adding up correctly')
end

% timescale for pressure decrease during aortic valve closure 
C = -decay_time / (R_distal * log(P_min/P_max)); 

fprintf("R_proximal = %.14f\n", R_proximal); 
fprintf("C = %.14f\n", C); 
fprintf("R_distal = %.14f\n", R_distal); 
R_total = R_proximal + R_distal
tau = C * R_distal


fprintf("right_pa_R_proximal = %.14f\n", R_proximal); 
fprintf("right_pa_C = %.14f\n", C); 
fprintf("right_pa_R_distal = %.14f\n", R_distal);
fprintf("left_pa_R_proximal = %.14f\n", R_proximal); 
fprintf("left_pa_C = %.14f\n", C); 
fprintf("left_pa_R_distal = %.14f\n", R_distal);


fig = figure; 

times_two_cycles = [t,t + cycle_length]; 
q_rv_exp = [flows_spline, flows_spline]; 

times = 0:dt:(2*cycle_length);
p_rv_exp  = Series_right_ventricle(times); 
p_pa_exp  = Series_pa(times); 

title('experimental pressures and flows')
subplot(2,1,1)
plot(times, p_rv_exp , 'k')
hold on
plot(times, p_pa_exp , ':k')
legend('P_{RV}', 'P_{PA}', 'Location','NorthEastOutside');
xlabel('t (s)')
ylabel('P (mmHg)')
subplot(2,1,2)
plot(times_two_cycles, q_rv_exp, 'k')
hold on
plot(times_two_cycles, zeros(size(q_rv_exp)), ':k')
axis([0 1.8 -150 250])
legend('Q RV', 'Location', 'NorthEastOutside')
xlabel('t (s)')
ylabel('Flow (ml/s), Net Flow (ml)')
set(fig, 'Position', [100, 100, 1000, 750])
set(fig,'PaperPositionMode','auto')
printfig(fig, 'bc_variables_experimental')

times_exp = times; 
save 'bc_variables_experimental.mat' times_two_cycles q_rv_exp times_exp p_rv_exp p_pa_exp









