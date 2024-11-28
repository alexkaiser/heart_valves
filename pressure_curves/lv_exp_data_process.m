
% quadrature spacing 
debug = true; 
if debug 
    dt = 5e-4; 
else 
    dt = 5e-6; 
end 

MMHG_TO_CGS = 1333.22368;

cycle_duration = 0.8; 

Q_goal_L_per_min = 5.6; 
Q_goal_ml_per_s = Q_goal_L_per_min * 1e3 / 60; 
Q_goal_ml_per_cycle = Q_goal_ml_per_s * cycle_duration; 

ejection_fraction_goal = 0.70;
end_diastolic_volume = Q_goal_ml_per_cycle / ejection_fraction_goal
ventricular_volume_initial = end_diastolic_volume - Q_goal_ml_per_cycle

% Poiseuille flow estimate on vessel 
L = 7; % cm 
radius = 1.25; % cm 
mu = 0.04; 
resistance_vessel = 1 / (8 * mu * L / (pi * radius^4))

% resistance mmHg 
resistance_lvot_mmHg = 0.0043; 
resistance_lvot = MMHG_TO_CGS * resistance_lvot_mmHg 

% estimate R ao via flow field differences at approximate peak 
% quite coarse but should not matter 
delta_p_av_mmHg = 7; 
q_av_est = 450; % ml/s
delta_p_av_dynescm2 = delta_p_av_mmHg * MMHG_TO_CGS; 

r_av = delta_p_av_dynescm2 / q_av_est

table = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_lv_la_ao.csv'); 

times_pressures = table.x; 
pressures_lv_raw = table.pressure_lv; 
pressures_aorta_raw = table.pressure_aorta; 
pressures_la_raw = table.pressure_la; 

pressure_raw_plot = false; 
if pressure_raw_plot
    figure; 
    plot(times_pressures, pressures_lv_raw)
    hold on 
    plot(times_pressures, pressures_la_raw)
    plot(times_pressures, pressures_aorta_raw)
end 

% normalize to cardiac cycle 
% start at zero 
times_pressures_adjusted = times_pressures(:) - times_pressures(1); 

% normalize times to desired cardiac cycle duration, both endpoints are included here 
times_pressures_adjusted = times_pressures_adjusted * cycle_duration / times_pressures_adjusted(end); 


bump_radius = 0.01; % 0.025; 
n_fourier_coeffs = 600; 
plots = false; 
t_mesh_one_cycle = 0:dt:cycle_duration; 
t = 0:dt:cycle_duration; 
% t = 0:dt:(cycle_duration*2); 

[a_0_pressure_lv, a_n_pressure_lv, b_n_pressure_lv, Series_pressure_lv, ~, ~, Series_pressure_lv_derivative] = ... 
    series_and_smooth([times_pressures_adjusted, pressures_lv_raw], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_pressure_lv = Series_pressure_lv(t); 
vals_series_pressure_lv_derivative = Series_pressure_lv_derivative(t); 

check_derivatives = false; 
if check_derivatives
    
    fd_derivatives = (Series_pressure_lv(t + dt) - Series_pressure_lv(t - dt)) / (2*dt); 
    
    figure; 
    plot(t, vals_series_pressure_lv_derivative)    
    hold on 
    plot(t, fd_derivatives)
    
    legend('series', 'finite diff derivative')
    
end 

% activation pressure proportional to ventricular pressure 
p_lv_activation_threshold = 20; 
activation_data_unscaled = (pressures_lv_raw > p_lv_activation_threshold) .* pressures_lv_raw;
activation_data = activation_data_unscaled / max(activation_data_unscaled);

[a_0_activation, a_n_activation, b_n_activation, Series_activation] = ... 
    series_and_smooth([times_pressures_adjusted, activation_data], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_activation = Series_activation(t);

base_name = 'fourier_coeffs';
suffix = '_lv_activation'
file_name = strcat(base_name, suffix, '.txt'); 
output_series_coeffs_to_txt(a_0_activation, a_n_activation, b_n_activation, n_fourier_coeffs, cycle_duration, file_name); 

[a_0_pressure_aorta, a_n_pressure_aorta, b_n_pressure_aorta, Series_pressure_aorta, ~, ~, Series_pressure_aorta_derivative] = ...
    series_and_smooth([times_pressures_adjusted, pressures_aorta_raw], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_pressure_aorta = Series_pressure_aorta(t); 
vals_series_pressure_aorta_derivative = Series_pressure_aorta_derivative(t); 


table_q_mitral = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_mitral.csv'); 
times_q_mitral = table_q_mitral.x; 
q_mitral_raw = table_q_mitral.q_mitral; 

times_q_mitral_adjusted = times_q_mitral - times_q_mitral(1);
times_q_mitral_adjusted = times_q_mitral_adjusted * cycle_duration / times_q_mitral_adjusted(end);

[a_0_q_mitral, a_n_q_mitral, b_n_q_mitral, Series_q_mitral, ~, ~, Series_q_mitral_derivative] = ...
    series_and_smooth([times_q_mitral_adjusted, q_mitral_raw], dt, bump_radius, n_fourier_coeffs, plots); 


vals_series_q_mitral = Series_q_mitral(t);
vals_series_q_mitral_derivative = Series_q_mitral_derivative(t);

table_q_aorta = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_aorta.csv'); 
times_q_aorta = table_q_aorta.x; 
q_aorta_raw = table_q_aorta.q_aorta; 

times_q_aorta_adjusted = times_q_aorta - times_q_aorta(1);
times_q_aorta_adjusted = times_q_aorta_adjusted * cycle_duration / times_q_aorta_adjusted(end);

[a_0_q_aorta, a_n_q_aorta, b_n_q_aorta, Series_q_aorta, ~, ~, Series_q_aorta_derivative] = ...
    series_and_smooth([times_q_aorta_adjusted, q_aorta_raw], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_q_aorta = Series_q_aorta(t);
vals_series_q_aorta_derivative = Series_q_aorta_derivative(t);


% normalize flow rate integrals 
% to be equal and to desired target
q_mitral_cumulative = dt * trapz(vals_series_q_mitral)
q_aorta_cumulative = dt * trapz(vals_series_q_aorta)

scaling_q_mitral = Q_goal_ml_per_cycle / q_mitral_cumulative;
scaling_q_aorta = Q_goal_ml_per_cycle / q_aorta_cumulative;

Series_q_mitral_scaled = @(t) scaling_q_mitral * Series_q_mitral(t);
Series_q_mitral_derivative_scaled = @(t) scaling_q_mitral * Series_q_mitral_derivative(t);

Series_q_aorta_scaled = @(t) scaling_q_aorta * Series_q_aorta(t);
Series_q_aorta_derivative_scaled = @(t) scaling_q_aorta * Series_q_aorta_derivative(t);

% re-compute tables with scaled series 
vals_series_q_mitral_scaled = Series_q_mitral_scaled(t);
vals_series_q_mitral_derivative_scaled = Series_q_mitral_derivative_scaled(t);

vals_series_q_aorta_scaled = Series_q_aorta_scaled(t);
vals_series_q_aorta_derivative_scaled = Series_q_aorta_derivative_scaled(t);

q_mitral_cumulative_scaled = dt * trapz(vals_series_q_mitral_scaled)
q_aorta_cumulative_scaled = dt * trapz(vals_series_q_aorta_scaled)


base_name = 'fourier_coeffs';
suffix = '_Q_mi'
file_name = strcat(base_name, suffix, '.txt'); 
output_series_coeffs_to_txt(scaling_q_mitral * a_0_q_mitral, scaling_q_mitral * a_n_q_mitral, scaling_q_mitral * b_n_q_mitral, n_fourier_coeffs, cycle_duration, file_name); 



% cgs units for pressure
Series_pressure_lv_cgs = @(t) MMHG_TO_CGS * Series_pressure_lv(t);
Series_pressure_lv_derivative_cgs = @(t) MMHG_TO_CGS * Series_pressure_lv_derivative(t);
vals_series_pressure_lv_cgs = Series_pressure_lv_cgs(t); 
vals_series_pressure_lv_derivative_cgs = Series_pressure_lv_derivative_cgs(t); 

Series_pressure_aorta_cgs = @(t) MMHG_TO_CGS * Series_pressure_aorta(t);
Series_pressure_aorta_derivative_cgs = @(t) MMHG_TO_CGS * Series_pressure_aorta_derivative(t);
vals_series_pressure_aorta_cgs = Series_pressure_aorta_cgs(t); 
vals_series_pressure_aorta_derivative_cgs = Series_pressure_aorta_derivative_cgs(t); 


vals_ventricular_volume = zeros(size(vals_series_q_aorta)); 
vals_ventricular_volume_deriv = zeros(size(vals_series_q_aorta)); 
vals_ventricular_volume(1) = ventricular_volume_initial;

for j=2:length(vals_ventricular_volume)
    vals_ventricular_volume(j) = vals_ventricular_volume(j-1) + dt * (vals_series_q_mitral_scaled(j) - vals_series_q_aorta_scaled(j)); 
end 

for j=1:length(vals_ventricular_volume)
    vals_ventricular_volume_deriv(j) = vals_series_q_mitral_scaled(j) - vals_series_q_aorta_scaled(j); 
end 

series_plots = false; 

if series_plots
    figure; 
    subplot(4,1,1)
    plot(t, vals_series_pressure_lv)
    hold on 
    plot(t, vals_series_pressure_aorta)

    subplot(4,1,2)
    hold on 
    plot(t, vals_series_q_mitral_scaled)
    plot(t, vals_series_q_aorta_scaled)

    subplot(4,1,3)
    plot(t, vals_ventricular_volume)

    subplot(4,1,4);
    plot(t, vals_series_activation);

    figure; 
    plot(t, vals_series_pressure_lv_derivative)
    hold on 
    plot(t, vals_series_pressure_aorta_derivative)
    plot(t, vals_series_q_mitral_derivative_scaled)
    plot(t, vals_series_q_aorta_derivative_scaled)
    plot(t, vals_ventricular_volume_deriv)
    title('derivatives')
    legend('p lv derivative', 'p aorta derivative', 'q mi derivative', 'q_ao_derivative', 'lv derivative')
end 



% run the lpn 
dt_lpn = 1e-4; 
n_cycles = 3; 
t_final = cycle_duration * n_cycles;
P_ao_initial = 94*MMHG_TO_CGS;
R_proximal = 83.6698220729; 
C =  0.00167055364456;
R_distal = 1287.64596307;

v_initial = ventricular_volume_initial;


KERCKHOFFS = true; 
Regazzoni_aaron = false; 

if KERCKHOFFS 
    Vrd = 26.1; % ml 
    Vrs = 18; % ml

    % parameters from KERCKHOFFS AMBE 2006 
    % lv 
    % note that min compliance is actually the larger value under confusing naming convention 
    C_min_ml_over_kPa = 11.0; 
    C_max_ml_over_kPa = 0.946; 

    % E_min_kPa = 1/C_min_ml_over_kPa; 
    % E_max_kPa = 1/C_max_ml_over_kPa; 

    ML_OVER_KPA_TO_ML_OVER_DYNEPERCM2 = 1e-4; 

    C_min_scaling = 5; 
    C_max_scaling = 5; 

    C_min_ml_over_dynespercm2 = C_min_scaling * C_min_ml_over_kPa * ML_OVER_KPA_TO_ML_OVER_DYNEPERCM2; 
    C_max_ml_over_dynespercm2 = C_max_scaling * C_max_ml_over_kPa * ML_OVER_KPA_TO_ML_OVER_DYNEPERCM2; 

    Emax = 1/C_max_ml_over_dynespercm2
    Emin = 1/C_min_ml_over_dynespercm2

    Emin_mmHg_over_ml = Emin / MMHG_TO_CGS 
    Emax_mmHg_over_ml = Emax / MMHG_TO_CGS

elseif Regazzoni_aaron
    
else 
    error('not implemented');
end 

    
    
R_av_closed = 100000;
steepness_av = 0.00001; 


[times_lpn, P_lv, Q_ao, P_ao, V_lv, R_tanh] = solve_lv_ao_lpn(dt_lpn, t_final, v_initial, Vrd, Vrs, Emax, Emin, ...
                                     Series_q_mitral_scaled, Series_activation, ...
                                     P_ao_initial, R_proximal, C, R_distal, r_av, R_av_closed, steepness_av);


figure; 
subplot(5,1,1)
plot(t, vals_series_pressure_lv)
hold on 
plot(t, vals_series_pressure_aorta)
plot(times_lpn, P_lv/MMHG_TO_CGS); 
plot(times_lpn, P_ao/MMHG_TO_CGS);
xlim([0 max(times_lpn)])

legend('lv exp', 'ao exp', 'lv lpn', 'ao lpn');


subplot(5,1,2)
hold on 
plot(t, vals_series_q_mitral_scaled)
plot(t, vals_series_q_aorta_scaled)
plot(times_lpn, Q_ao)
xlim([0 max(times_lpn)])

legend('q mi exp', 'q ao exp', 'q ao lpn')

subplot(5,1,3)
plot(t, vals_ventricular_volume)
plot(times_lpn, V_lv)
xlim([0 max(times_lpn)])
legend('V integrated', 'V lpn')

subplot(5,1,4);
plot(t, vals_series_activation);
xlim([0 max(times_lpn)])

subplot(5,1,5)
plot(times_lpn, R_tanh)
xlim([0 max(times_lpn)])

figure; 
hold on 
plot(V_lv, P_lv/MMHG_TO_CGS)
xlabel('V ml')
ylabel('P LV mmHg')










output_to_sv0d = false; 
if output_to_sv0d
    % output 
    format long 

    f = fopen("array_values.txt", "w");

    % fprintf("time\n");
    % for j = 1:length(t)
    %     fprintf(f, "%.14f, ", t(j));
    % end 
    % fprintf(f, "\n");
    % 
    % fprintf(f, "vals_series_q_aorta_scaled\n");
    % for j = 1:length(t)
    %     fprintf(f, "%.14f, ", vals_series_q_aorta_scaled(j));
    % end 
    % fprintf(f, "\n");
    % 
    % fprintf(f, "vals_series_q_mitral_scaled\n");
    % for j = 1:length(t)
    %     fprintf(f, "%.14f, ", vals_series_q_mitral_scaled(j));
    % end 
    % fprintf(f, "\n");
    % 
    % fprintf(f, "vals_series_pressure_lv_cgs\n")
    % for j = 1:length(t)
    %     fprintf(f, "%.14f, ", vals_series_pressure_lv_cgs(j));
    % end 
    % fprintf(f, "\n");
    % 
    % fprintf(f, "vals_series_pressure_aorta_cgs\n")
    % for j = 1:length(t)
    %     fprintf(f, "%.14f, ", vals_series_pressure_aorta_cgs(j));
    % end 
    % fprintf(f, "\n");


    % 
    fprintf(f, '        "y": {\n');
    % fprintf(f, '        "flow:ventricle:valve1": [\n' );
    % fprintf(f, '            ');
    % for j = 1:length(t)
    %     fprintf(f, '%.14f', vals_series_q_aorta_scaled(j));
    %     if j < length(t)
    %         fprintf(f, ', ');
    %     end 
    % end 
    % fprintf(f, '\n');
    % fprintf(f, '        ],');

    print_var_string(f,t,'flow:ventricle:valve1', vals_series_q_aorta_scaled)
    print_var_string(f,t,'pressure:ventricle:valve1', vals_series_pressure_lv_cgs)
    print_var_string(f,t,'flow:vessel:OUTLET', vals_series_q_aorta_scaled)
    print_var_string(f,t,'pressure:vessel:OUTLET', vals_series_pressure_aorta_cgs)
    print_var_string(f,t,'Vc:ventricle', vals_ventricular_volume);

    fprintf(f, '    },\n');

    fprintf(f, '    "dy": {\n');
    print_var_string(f,t,'flow:ventricle:valve1', vals_series_q_aorta_derivative_scaled)
    print_var_string(f,t,'pressure:ventricle:valve1', vals_series_pressure_lv_derivative_cgs)
    print_var_string(f,t,'flow:vessel:OUTLET', vals_series_q_aorta_derivative_scaled)
    print_var_string(f,t,'pressure:vessel:OUTLET', vals_series_pressure_aorta_derivative_cgs)
    print_var_string(f,t,'Vc:ventricle', vals_ventricular_volume_deriv);
    fprintf(f, '    },\n');

end 

function print_var_string(f,t,name,vals)

    fprintf(f, '        "%s": [\n', name);
    fprintf(f, '            ');
    for j = 1:length(t)
        fprintf(f, '%.14f', vals(j));
        if j < length(t)
            fprintf(f, ', ');
        end 
    end 
    fprintf(f, '\n');
    fprintf(f, '        ],\n');

end 







