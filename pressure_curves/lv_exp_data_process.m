
% quadrature spacing 
debug = false; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 

MMHG_TO_CGS = 1333.22368;




historical_3 = true;
if historical_3
    
    HR = 85;
    
    cycle_duration = 60/HR
    
    % VTI version 
    % SV = 71.70;
    
    LVEDV = 121.7 
    LVESV = 56.9

    % SV = LVEDV - LVESV
    SV = 60

    % from echo 
    EF_echo = .58 
    EF_computed = SV / LVEDV

    % use clinical value 
    Q_goal_ml_per_cycle = SV; 

    ventricular_volume_initial = LVESV;
    
    adjust_venctricular_volume_by_start = true; 
    start_time_in_cycle = 0.1;
    
    % scales the waveform to these values 
    p_systolic_scaling = 105;
    p_diastolic_scaling = 65;
    
    basename_suffix = '_hist_3';
    
else    
    % basic case  
    cycle_duration = 0.8;
    Q_goal_L_per_min = 5.6; 
    ejection_fraction_goal = 0.70;
    
    Q_goal_ml_per_s = Q_goal_L_per_min * 1e3 / 60; 
    Q_goal_ml_per_cycle = Q_goal_ml_per_s * cycle_duration; 


    end_diastolic_volume = Q_goal_ml_per_cycle / ejection_fraction_goal
    ventricular_volume_initial = end_diastolic_volume - Q_goal_ml_per_cycle
    
    % permutes values to start at this time 
    % time arrays still start at zero  
    start_time_in_cycle = 0.15;
    
    basename_suffix = '';
end 
    
base_name = strcat('fourier_coeffs', basename_suffix);




% Poiseuille flow estimate on vessel 
L = 7; % cm 
radius = 1.25; % cm 
mu = 0.04; 
resistance_vessel = 1 / (8 * mu * L / (pi * radius^4));

% resistance mmHg 
resistance_lvot_mmHg = 0.0043; 
resistance_lvot = MMHG_TO_CGS * resistance_lvot_mmHg;

% estimate R ao via flow field differences at approximate peak 
% quite coarse but should not matter 
delta_p_av_mmHg = 7; 
q_av_est = 450; % ml/s
delta_p_av_dynescm2 = delta_p_av_mmHg * MMHG_TO_CGS; 

r_av = delta_p_av_dynescm2 / q_av_est;


two_hill = true; 
if two_hill
    % hill function parameters 
    % brown AMBE 2023 
%     tau_1 = 0.0725 * 2 * cycle_duration 
%     tau_2 = 0.4503 * cycle_duration
% 
%     m1 = 2.7463
%     m2 = 21.5683

    if historical_3
        t_shift =  0.2592964132032478;
        tau_1 =  0.13421838866217944;
        tau_2 =  0.327400222418779;
        m1 =  1.6204463203333281;
        m2 =  18.37110863697137;
    else 
        t_shift =  0.4080694895094201;
        tau_1 =  0.7999999995665112;
        tau_2 =  0.37189471760890913;
        m1 =  1.0916763390925226;
        m2 =  18.491969573464523;

    end 



%     g1 = @(t) (t > 0) .* (t./tau_1).^m1; 
%     g2 = @(t) (t > 0) .* (t./tau_2).^m2; 


    g1 = @(t) (t./tau_1).^m1; 
    g2 = @(t) (t./tau_2).^m2; 
    
    r1 = @(t) g1(t) ./ (1 + g1(t)); 
    r2 = @(t) 1 ./ (1 + g2(t)); 
    two_hill_product = @(t) (g1(t) ./ (1 + g1(t))) .* (1 ./ (1 + g2(t))); 

    % broken 
    % two_hill_product_maximum_opt = -fminbnd(@(t) -two_hill_product(t), 0, cycle_duration); 

    % hack opt 
    times_hill = 0:dt:cycle_duration; % linspace(0,cycle_duration,cycle_duration/dt); 
    two_hill_product_maximum = max(two_hill_product(times_hill)); 

    % times_hill = linspace(0,cycle_duration,cycle_duration/dt); 
    
    k_coeff_two_hill = 1 / two_hill_product_maximum; 

    t_shift_hill = t_shift;
    
    % no periodicity handled here 
    % two_hill_function = @(t) k_coeff_two_hill * two_hill_product(t - t_shift_hill); 

    two_hill_function = @(t) k_coeff_two_hill * two_hill_product(mod(t - t_shift_hill,cycle_duration)); 
    
    % two_hill_function = @(t) two_hill_product(t); 

    figure; 
%     plot(times_hill, two_hill_product(times_hill)); 
    hold on    
    plot(times_hill, two_hill_function(times_hill)); 
    
%     plot(times_hill, g1(times_hill))
%     plot(times_hill, max(g1(times_hill), zeros(size(times_hill))) )
%     plot(times_hill, g1(times_hill)/max(abs(g1(times_hill))))
%     plot(times_hill, g2(times_hill)/max(abs(g2(times_hill))))

    plot(times_hill, r1(times_hill))
    plot(times_hill, r2(times_hill))

    % plot(times_hill, two_hill_function_periodic(times_hill)); 
    
    legend('hill function', 'r1', 'r2')
    
    
    title('hill function')
    
end 




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
    series_and_smooth([times_pressures_adjusted, pressures_lv_raw], dt, bump_radius, n_fourier_coeffs, plots, start_time_in_cycle); 

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

% no time shift here, two hill includes the time shift from optimizer
[a_0_activation_two_hill, a_n_activation_two_hill, b_n_activation_two_hill, Series_activation_two_hill] = ... 
    series_and_smooth([t', two_hill_function(t)'], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_activation_two_hill = Series_activation_two_hill(t);

figure; 
plot(t, vals_series_activation_two_hill)
hold on 
plot(t', two_hill_function(t)', 'k--')
legend('series', 'inputs')

title('two hill act from series in lv exp data process')

suffix = '_lv_activation_two_hill';
file_name = strcat(base_name, suffix, '.txt'); 
output_series_coeffs_to_txt(a_0_activation_two_hill, a_n_activation_two_hill, b_n_activation_two_hill, n_fourier_coeffs, cycle_duration, file_name); 



% activation pressure proportional to ventricular pressure 
% turn this off 
% p_lv_activation_threshold = 20; 
% activation_data_unscaled = (pressures_lv_raw > p_lv_activation_threshold) .* pressures_lv_raw;
% activation_data = activation_data_unscaled / max(activation_data_unscaled);
% 
% [a_0_activation, a_n_activation, b_n_activation, Series_activation] = ... 
%     series_and_smooth([times_pressures_adjusted, activation_data], dt, bump_radius, n_fourier_coeffs, plots, start_time_in_cycle); 
% 
% vals_series_activation = Series_activation(t);
% 
% 
% suffix = '_lv_activation'
% file_name = strcat(base_name, suffix, '.txt'); 
% output_series_coeffs_to_txt(a_0_activation, a_n_activation, b_n_activation, n_fourier_coeffs, cycle_duration, file_name); 

[a_0_pressure_aorta, a_n_pressure_aorta, b_n_pressure_aorta, Series_pressure_aorta, ~, ~, Series_pressure_aorta_derivative] = ...
    series_and_smooth([times_pressures_adjusted, pressures_aorta_raw], dt, bump_radius, n_fourier_coeffs, plots, start_time_in_cycle); 

vals_series_pressure_aorta = Series_pressure_aorta(t); 
vals_series_pressure_aorta_derivative = Series_pressure_aorta_derivative(t); 

if exist('p_systolic_scaling', 'var') && exist('p_diastolic_scaling', 'var') 
    
    p_sys_orig = max(vals_series_pressure_aorta);
    p_dia_orig = min(vals_series_pressure_aorta);
    
    pulse_pressure_orig = p_sys_orig - p_dia_orig;
    
    pulse_pressure_goal = p_systolic_scaling - p_diastolic_scaling;
    
    pressure_scaling_ao = pulse_pressure_goal / pulse_pressure_orig;
    
    pressure_shift_mean = pressure_scaling_ao * p_dia_orig - p_diastolic_scaling;
    
    a_0_pressure_aorta_adjust = pressure_scaling_ao * a_0_pressure_aorta - pressure_shift_mean; 
    a_n_pressure_aorta_adjust = pressure_scaling_ao * a_n_pressure_aorta;
    b_n_pressure_aorta_adjust = pressure_scaling_ao * b_n_pressure_aorta;
        
    aorta_series_no_array = @(t) a_0_pressure_aorta_adjust + sum(a_n_pressure_aorta_adjust .* cos((2*pi/cycle_duration) * (1:n_fourier_coeffs) .* t)' + ...  
                                 b_n_pressure_aorta_adjust .* sin((2*pi/cycle_duration) * (1:n_fourier_coeffs) .* t)' );   

    Series_pressure_aorta_adjust = @(t) arrayfun(aorta_series_no_array, t); 
    
    vals_series_pressure_aorta_adjust = Series_pressure_aorta_adjust(t);

    p_sys_adjust = max(vals_series_pressure_aorta_adjust);
    p_dia_adjust = min(vals_series_pressure_aorta_adjust);
    
    pulse_pressure_adjust = p_sys_adjust - p_dia_adjust;
    
    % scale by multiple of aortic systolic pressure 
    pressure_scaling_lv = p_sys_adjust / p_sys_orig;
    
    a_0_pressure_lv_adjust = pressure_scaling_lv * a_0_pressure_lv; 
    a_n_pressure_lv_adjust = pressure_scaling_lv * a_n_pressure_lv;
    b_n_pressure_lv_adjust = pressure_scaling_lv * b_n_pressure_lv;
        
    lv_series_no_array = @(t) a_0_pressure_lv_adjust + sum(a_n_pressure_lv_adjust .* cos((2*pi/cycle_duration) * (1:n_fourier_coeffs) .* t)' + ...  
                                 b_n_pressure_lv_adjust .* sin((2*pi/cycle_duration) * (1:n_fourier_coeffs) .* t)' );   

    Series_pressure_lv_adjust = @(t) arrayfun(lv_series_no_array, t); 
    
    vals_series_pressure_lv_adjust = Series_pressure_lv_adjust(t);
    
    figure; 
    plot(t, vals_series_pressure_aorta)
    hold on 
    plot(t, vals_series_pressure_aorta_adjust)
    plot(t, vals_series_pressure_lv)
    plot(t, vals_series_pressure_lv_adjust)
        
    p_lv_max = max(vals_series_pressure_lv);
    p_lv_min = min(vals_series_pressure_lv);
    p_lv_max_adjust = max(vals_series_pressure_lv_adjust);
    p_lv_min_adjust = min(vals_series_pressure_lv_adjust);

    % just reset the arrays 
    a_0_pressure_aorta = a_0_pressure_aorta_adjust;
    a_n_pressure_aorta = a_n_pressure_aorta_adjust;
    b_n_pressure_aorta = b_n_pressure_aorta_adjust;
    Series_pressure_aorta = Series_pressure_aorta_adjust;
    
    a_0_pressure_lv = a_0_pressure_lv_adjust;
    a_n_pressure_lv = a_n_pressure_lv_adjust;
    b_n_pressure_lv = b_n_pressure_lv_adjust;
    Series_pressure_lv = Series_pressure_lv_adjust;
    
    vals_series_pressure_aorta = vals_series_pressure_aorta_adjust;
    vals_series_pressure_lv = vals_series_pressure_lv_adjust;


end 



table_q_mitral = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_mitral.csv'); 
times_q_mitral = table_q_mitral.x; 
q_mitral_raw = table_q_mitral.q_mitral; 

times_q_mitral_adjusted = times_q_mitral - times_q_mitral(1);
times_q_mitral_adjusted = times_q_mitral_adjusted * cycle_duration / times_q_mitral_adjusted(end);

[a_0_q_mitral, a_n_q_mitral, b_n_q_mitral, Series_q_mitral, ~, ~, Series_q_mitral_derivative] = ...
    series_and_smooth([times_q_mitral_adjusted, q_mitral_raw], dt, bump_radius, n_fourier_coeffs, plots, start_time_in_cycle); 


vals_series_q_mitral = Series_q_mitral(t);
vals_series_q_mitral_derivative = Series_q_mitral_derivative(t);

table_q_aorta = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_aorta.csv'); 
times_q_aorta = table_q_aorta.x; 
q_aorta_raw = table_q_aorta.q_aorta; 

times_q_aorta_adjusted = times_q_aorta - times_q_aorta(1);
times_q_aorta_adjusted = times_q_aorta_adjusted * cycle_duration / times_q_aorta_adjusted(end);

[a_0_q_aorta, a_n_q_aorta, b_n_q_aorta, Series_q_aorta, ~, ~, Series_q_aorta_derivative] = ...
    series_and_smooth([times_q_aorta_adjusted, q_aorta_raw], dt, bump_radius, n_fourier_coeffs, plots, start_time_in_cycle); 

vals_series_q_aorta = Series_q_aorta(t);
vals_series_q_aorta_derivative = Series_q_aorta_derivative(t);


% normalize flow rate integrals 
% to be equal and to desired target
q_mitral_cumulative = dt * trapz(vals_series_q_mitral);
q_aorta_cumulative = dt * trapz(vals_series_q_aorta);

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



suffix = '_Q_mi';
file_name = strcat(base_name, suffix, '.txt'); 
output_series_coeffs_to_txt(scaling_q_mitral * a_0_q_mitral, scaling_q_mitral * a_n_q_mitral, scaling_q_mitral * b_n_q_mitral, n_fourier_coeffs, cycle_duration, file_name); 


if exist('adjust_venctricular_volume_by_start', 'var') && adjust_venctricular_volume_by_start
    % adds fraction of Q_mi that occurs during the shift to the initial ventricular volume 

    t_final_adjust = t(end);
    t_start_adjust = t_final_adjust - start_time_in_cycle;
    min_idx = find(t >= t_start_adjust,1);
    adjust_volume = dt * trapz(vals_series_q_mitral_scaled(min_idx:end))
    ventricular_volume_initial = ventricular_volume_initial + adjust_volume
end 



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

series_plots = true; 

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
    plot(t, vals_series_activation_two_hill);

    %     figure; 
    %     plot(t, vals_series_pressure_lv_derivative)
    %     hold on 
    %     plot(t, vals_series_pressure_aorta_derivative)
    %     plot(t, vals_series_q_mitral_derivative_scaled)
    %     plot(t, vals_series_q_aorta_derivative_scaled)
    %     plot(t, vals_ventricular_volume_deriv)
    %     title('derivatives')
    %     legend('p lv derivative', 'p aorta derivative', 'q mi derivative', 'q_ao_derivative', 'lv derivative')
end 


run_matlab_zerod = false; 
if run_matlab_zerod
    % run the lpn 
    dt_lpn = 1e-4; 
    n_cycles = 1; 
    t_final = cycle_duration * n_cycles;
    P_ao_initial = 94*MMHG_TO_CGS;
    R_proximal = 83.6698220729; 
    C =  0.00167055364456;
    R_distal = 1287.64596307;

    v_initial = ventricular_volume_initial;


    KERCKHOFFS = true; 

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

    else 
        error('not implemented');
    end 



    R_av_closed = 100000;
    steepness_av = 0.00001; 


    [times_lpn, P_lv, Q_ao, P_ao, V_lv, R_tanh] = solve_lv_ao_lpn(dt_lpn, t_final, v_initial, Vrd, Vrs, Emax, Emin, ...
                                         Series_q_mitral_scaled, Series_activation, ...
                                         P_ao_initial, R_proximal, C, R_distal, r_av, R_av_closed, steepness_av);
 


    figure;
    subplot(6,1,1)
    plot(t, vals_series_pressure_lv)
    hold on 
    plot(t, vals_series_pressure_aorta)
    plot(times_lpn, P_lv/MMHG_TO_CGS); 
    plot(times_lpn, P_ao/MMHG_TO_CGS);
    xlim([0 max(times_lpn)])

    legend('lv exp', 'ao exp', 'lv lpn', 'ao lpn');


    subplot(6,1,2)
    hold on 
    plot(t, vals_series_q_mitral_scaled)
    plot(t, vals_series_q_aorta_scaled)
    plot(times_lpn, Q_ao)
    xlim([0 max(times_lpn)])
    legend('q mi exp', 'q ao exp', 'q ao lpn')


    subplot(6,1,3)
    plot(times_hill, two_hill_function(times_hill)); 


    subplot(6,1,4)
    plot(t, vals_ventricular_volume)
    hold on 
    plot(times_lpn, V_lv)
    xlim([0 max(times_lpn)])
    legend('V integrated', 'V lpn')

    subplot(6,1,5);
    plot(t, vals_series_activation_two_hill);
    hold on 
    plot(t, vals_series_activation);
    xlim([0 max(times_lpn)])
    legend('two hill', 'act from p lv')
    title('act')

    subplot(6,1,6)
    plot(times_lpn, R_tanh)
    xlim([0 max(times_lpn)])
    title('r valve')
    
    figure; 
    hold on 
    plot(V_lv, P_lv/MMHG_TO_CGS)
    xlabel('V ml')
    ylabel('P LV mmHg')

end


p_ao_mean_mmHg = mean(vals_series_pressure_aorta)
p_ao_mean_cgs  = mean(vals_series_pressure_aorta_cgs)

p_ao_initial_mmHg = vals_series_pressure_aorta(1)
p_ao_initial_cgs  = vals_series_pressure_aorta_cgs(1)




% % estimate linear pressure volume relation 
% t_0 = 0; 
% t_1 = .2; 
% 
% t_0_idx = find(times >= t_0, 1);
% t_1_idx = find(times >= t_1, 1);
% 
% vol_0 = vals_ser






output_to_sv0d = true; 
if output_to_sv0d
    % output 
    format long 

    f = fopen("array_values.json", "w");
    
    % 
    fprintf(f,'"Vc:ventricle": %.14f,\n', vals_ventricular_volume(1));
    fprintf(f,'\n');
    
    % times
    print_var_string(f,t,'Q',vals_series_q_mitral_scaled);
    print_var_string(f,t,'t',t);
    fprintf(f,'\n');
    
   
    % 
    fprintf(f, '        "y": {\n');
    print_var_string(f,t,'flow:ventricle:valve1', vals_series_q_aorta_scaled)
    print_var_string(f,t,'pressure:ventricle:valve1', vals_series_pressure_lv_cgs)
    print_var_string(f,t,'flow:vessel:OUTLET', vals_series_q_aorta_scaled)
    print_var_string(f,t,'pressure:vessel:OUTLET', vals_series_pressure_aorta_cgs)
    print_var_string(f,t,'Vc:ventricle', vals_ventricular_volume);

    fprintf(f, '    },\n');

    plot_pressure_checks = true; 
    if plot_pressure_checks
        figure; 
        hold on; 
        plot(t, vals_series_pressure_lv)
        plot(t, vals_series_pressure_lv_cgs/MMHG_TO_CGS)

        plot(t, vals_series_pressure_aorta)
        plot(t, vals_series_pressure_aorta_cgs/MMHG_TO_CGS)

        legend('lv mmHg', 'lv converted', 'ao mmHg', 'ao converted')

    end 


%     fprintf(f, '    "dy": {\n');
%     print_var_string(f,t,'flow:ventricle:valve1', vals_series_q_aorta_derivative_scaled)
%     print_var_string(f,t,'pressure:ventricle:valve1', vals_series_pressure_lv_derivative_cgs)
%     print_var_string(f,t,'flow:vessel:OUTLET', vals_series_q_aorta_derivative_scaled)
%     print_var_string(f,t,'pressure:vessel:OUTLET', vals_series_pressure_aorta_derivative_cgs)
%     print_var_string(f,t,'Vc:ventricle', vals_ventricular_volume_deriv);
%     fprintf(f, '    },\n');

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







