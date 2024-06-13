

MMHG_TO_CGS = 1333.22368;

cycle_length = 0.8; 


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

Emin_mmHg_over_ml = Emin / MMHG_TO_CGS; 
Emax_mmHg_over_ml = Emax / MMHG_TO_CGS;


% estimate R ao 
delta_p_av_mmHg = 7; 
q_av_est = 450; % ml/s
delta_p_av_dynescm2 = delta_p_av_mmHg * MMHG_TO_CGS; 

r_av = delta_p_av_dynescm2 / q_av_est; 

cos_power = 1; 


Vrd = 26.1; % ml 
Vrs = 18.0; % ml
t_active = 0.6; % s 
t_twitch = 0.2; % s 

% 
inductance = 0.000351787; 



two_hill = false; 
if two_hill
    % hill function parameters 
    % brown AMBE 2023 
    tau_1 = 0.0725 * cycle_length; 
    tau_2 = 0.6 * 0.4503 * cycle_length; 
    m1 = 2.7463; 
    m2 = 21.5683; 

    g1 = @(t) (t./tau_1).^m1; 
    g2 = @(t) (t./tau_2).^m2; 

    two_hill_product = @(t) (g1(t) ./ (1 + g1(t))) .* (1 ./ (1 + g2(t))); 

    % broken 
    % two_hill_product_maximum_opt = -fminbnd(@(t) -two_hill_product(t), 0, cycle_length); 

    % hack opt 
    times_hill = linspace(0,cycle_length,1e5); 
    two_hill_product_maximum = max(two_hill_product(times_hill)); 


    k_coeff_two_hill = (Emax - Emin) / two_hill_product_maximum; 

    two_hill_function = @(t) k_coeff_two_hill * two_hill_product(t) + Emin; 

    % two_hill_function = @(t) two_hill_product(t); 

    figure; 
    times_hill = linspace(0,cycle_length,1e4); 
    plot(times_hill, two_hill_product(times_hill)); 
    title('hill function product')

    figure; 
    times_hill = linspace(0,cycle_length,1e4); 
    plot(times_hill, two_hill_function(times_hill)); 
    title('hill function')
end 


Q_goal_L_per_min = 5.6; 
Q_goal_ml_per_s = Q_goal_L_per_min * 1e3 / 60; 
Q_goal_ml_per_cycle = Q_goal_ml_per_s * cycle_length; 

% quadrature spacing 
debug = false; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 


points_one_cycle_Q_mi = [
-0.74695,-0.131
-0.72432,30.817
-0.712413,68.579
-0.699234,139.489
-0.69193,181.232
-0.68599,201.827
-0.676983,220.149
-0.667771,211.615
-0.653874,188.242
-0.639976,164.869
-0.603057,121.016
-0.557017,80.629
-0.51876,61.354
-0.459132,38.163
-0.416364,26.906
-0.382893,35.611
-0.360181,55.702
-0.352777,84.303
-0.348474,119.749
-0.342555,143.201
-0.329054,171.827
-0.318165,143.299
-0.313414,119.889
-0.307038,83.343
-0.302243,54.22
-0.297495,31.381
-0.292678,-0.599
-0.286337,-32.574
-0.277017,-55.394
-0.27382,-74.81
-0.267649,-84.499
-0.261621,-75.332
-0.255846,-33.023
-0.25457,-0.447
-0.254452,-15.875
-0.253163,14.987
-0.250241,31.57
-0.243995,12.167
-0.240855,0.179
-0.239291,-4.957
-0.23332,11.638
-0.225607,-0.331
-0.19512,-0.209
-0.137187,-1.12
-0.100607,-0.402
-0.079266,-0.317
-0.062498,-0.25]; 

points_one_cycle_p_lv_mmHg = [
-0.746864,29.758
-0.733269,17.538
-0.710503,4.179
-0.684491,2.751
-0.646184,4.161
-0.606373,3.865
-0.554301,4.135
-0.497627,4.971
-0.424072,8.075
-0.365836,10.899
-0.321378,14.012
-0.289215,14.287
-0.265991,30.757
-0.256468,52.63
-0.251674,65.696
-0.243734,84.16
-0.240415,100.921
-0.234028,117.964
-0.218636,123.073
% -0.191049,124.486
-0.160369,127.886
-0.132791,128.731
-0.114445,126.737
-0.102285,120.768
-0.097817,112.528
-0.0903,103.435
-0.08584,94.627
-0.081416,83.546
-0.075451,73.033
-0.068008,59.11
-0.063561,49.45
-0.057597,38.937
]; 


% start at zero 
times_Q_mi_adjusted = points_one_cycle_Q_mi(:,1) - points_one_cycle_Q_mi(1,1); 

% normalize times to desired cardiac cycle duration 
times_Q_mi_adjusted = times_Q_mi_adjusted * cycle_length / times_Q_mi_adjusted(end); 

Q_mi_adjusted = points_one_cycle_Q_mi(:,2);  

points_one_cycle_Q_mi_adjusted = [times_Q_mi_adjusted, Q_mi_adjusted]; 


% start at zero 
times_p_lv_mmHg_adjusted = points_one_cycle_p_lv_mmHg(:,1) - points_one_cycle_p_lv_mmHg(1,1); 

% normalize times to desired cardiac cycle duration 
times_p_lv_mmHg_adjusted = times_p_lv_mmHg_adjusted * cycle_length / times_p_lv_mmHg_adjusted(end); 

p_lv_mmHg_adjusted = points_one_cycle_p_lv_mmHg(:,2);  

% threshold 
p_lv_activation_threshold = 20; 
p_lv_mmHg_adjusted = (p_lv_mmHg_adjusted > p_lv_activation_threshold) .* p_lv_mmHg_adjusted; 
p_lv_normalized = p_lv_mmHg_adjusted / max(p_lv_mmHg_adjusted); 

points_one_cycle_lv_activation= [times_p_lv_mmHg_adjusted, p_lv_normalized]; 



% plot(times_Q_mi_adjusted, points_one_cycle_Q_mi(:,2)); 

base_name = 'fourier_coeffs';
suffix = '_Q_mi'

file_name = strcat(base_name, suffix, '.txt'); 

bump_radius = .025; 
n_fourier_coeffs = 600; 
plots = false; 

[a_0_Q_mi a_n_Q_mi b_n_Q_mi Series_Q_mi ] = series_and_smooth(points_one_cycle_Q_mi_adjusted, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_Q_mi_series = Series_Q_mi(t); 

% fig = figure; 
% plot(t, vals_Q_mi_series, 'k'); 
% hold on
% title('Q mi')
% xlabel('t')
% ylabel('Q (ml/s)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')

% trap integral 
Q_mi_total = sum(vals_Q_mi_series(1:(end-1)) * dt); 

scaling_q_mi = Q_goal_ml_per_cycle / Q_mi_total; 

% scale mitral to goal 
Series_Q_mi_scaled = @(t) scaling_q_mi * Series_Q_mi(t); 

output_series_coeffs_to_txt(scaling_q_mi * a_0_Q_mi, scaling_q_mi * a_n_Q_mi, scaling_q_mi * b_n_Q_mi, n_fourier_coeffs, cycle_length, file_name); 



base_name = 'fourier_coeffs';
suffix = '_lv_activation'

file_name = strcat(base_name, suffix, '.txt'); 

bump_radius = .025; 
n_fourier_coeffs = 600; 
plots = false; 

[a_0_lv_act, a_n_lv_act, b_n_lv_act, Series_lv_activation] = series_and_smooth(points_one_cycle_lv_activation, dt, bump_radius, n_fourier_coeffs, plots); 

output_series_coeffs_to_txt(a_0_lv_act, a_n_lv_act, b_n_lv_act, n_fourier_coeffs, cycle_length, file_name); 

t = 0:dt:cycle_length; 
vals_lv_activation = Series_lv_activation(t); 

fig = figure; 
plot(t, vals_lv_activation, 'k'); 
hold on
title('lv act')
xlabel('t')
ylabel('act')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')

% % trap integral 
% Q_mi_total = sum(vals_Q_mi_series(1:(end-1)) * dt); 
% 
% % scale mitral to goal 
% Series_Q_mi_scaled = @(t) (Q_goal_ml_per_cycle / Q_mi_total) * Series_Q_mi(t); 






load('bc_data_aortic_20978583_384_63b7a5b.mat', 'times', 'q_aorta', 'p_aorta'); 
times_ao = times; 

min_idx_times_ao = find(times >= cycle_length, 1); 
max_idx_times_ao = find(times > (2*cycle_length), 1) - 1; 

times_cycle_2   = times_ao(min_idx_times_ao:max_idx_times_ao); 
q_aorta_cycle_2 = q_aorta(min_idx_times_ao:max_idx_times_ao); 

dt_Q_ao_series = times_cycle_2(2) - times_cycle_2(1); 
[~, ~, ~, Series_Q_Ao] = fourier_series_uniform(times_cycle_2, q_aorta_cycle_2, cycle_length, n_fourier_coeffs, dt_Q_ao_series); 

vals_Q_Ao_series = Series_Q_Ao(t); 

Q_Ao_total = sum(vals_Q_Ao_series(1:(end-1)) * dt); 

% plot(t, vals_Q_Ao_series, 'k'); 

% scale Ao too for conservation 
Series_Q_Ao_scaled = @(t) (Q_goal_ml_per_cycle / Q_Ao_total) * Series_Q_Ao(t); 


p_aorta_cycle_2 = p_aorta(min_idx_times_ao:max_idx_times_ao); 

[~, ~, ~, Series_P_Ao_mmHg] = fourier_series_uniform(times_cycle_2, p_aorta_cycle_2, cycle_length, n_fourier_coeffs, dt_Q_ao_series); 

Series_P_Ao_dynescm2 = @(t) MMHG_TO_CGS * Series_P_Ao_mmHg(t); 

dt = 1e-3; 

t_offset_Q_mi = 0; 

n_cycles = 3; 
max_t = n_cycles * cycle_length; 

times = 0:dt:max_t; 

n_times = length(times); 

P_lv_in  = zeros(n_times, 1); 
P_lv_out = zeros(n_times, 1); 
V_lv     = zeros(n_times, 1); 
act      = zeros(n_times, 1); 
Vrest    = zeros(n_times, 1); 
Elas     = zeros(n_times, 1); 

% start at diastolic rest volume 
V_lv(1) = Vrd; 

% both pressures can start at zero 

% forcing terms 
Q_mi = Series_Q_mi_scaled(times - t_offset_Q_mi); 
act  = Series_lv_activation(times - t_offset_Q_mi); 

prescribe_Q_ao = false; 
if prescribe_Q_ao
    Q_ao = Series_Q_Ao_scaled(times); 
else 
    Q_ao = zeros(n_times, 1); 
    P_ao_distal = Series_P_Ao_dynescm2(times); 
end 

V_lv(1) = Vrd; 

for n = 1:(n_times-1)

    t = dt*n; 

    if ~prescribe_Q_ao
        if P_lv_out(n) > P_ao_distal(n)
            Q_ao(n) = (1/r_av) * (P_lv_out(n) - P_ao_distal(n)); 
        else 
            Q_ao(n) = 0; 
        end 
    end 
    
    % eqn 3 for ventricular volume 
    V_lv(n+1) = V_lv(n) + dt * (Q_mi(n) - Q_ao(n)); 

    % calculate current elastance Elas 
    t_in_cycle = mod(t, cycle_length);

    t_contract = 0;
    if (t_in_cycle >= t_active)
        t_contract = t_in_cycle - t_active;
    end

%     act(n+1) = 0;
%     if (t_contract <= t_twitch)
%         % act(n+1) = (-0.5 * cos(2 * pi * t_contract / t_twitch) + 0.5)^cos_power;
%     end 

    if t < 0.2
        act_temp = 0; 
    else 
        act_temp = act(n+1); 
    end 

    Vrest(n+1) = (1.0 - act_temp) * (Vrd - Vrs) + Vrs;
    Elas(n+1) = (Emax - Emin) * act_temp+ Emin;

    % eqn 2 
    % check the time lags here 
    P_lv_in(n+1) = Elas(n+1) * (V_lv(n+1) - Vrest(n+1)); 

    if n > 1
        P_lv_out(n+1) = P_lv_in(n+1) - inductance * (Q_ao(n) - Q_ao(n-1))/dt; 
    else 
        P_lv_out(n+1) = P_lv_in(n+1); 
    end 
    
end 






P_lv_out_mmHg = P_lv_out / MMHG_TO_CGS;
P_lv_in_mmHg  = P_lv_in  / MMHG_TO_CGS; 

fig = figure; 
subplot(5,1,1)
hold on 
plot(times, Q_mi)
plot(times, Q_ao)
ylabel('ml')
legend('Q mi', 'Q ao')


subplot(5,1,2)
hold on 
plot(times, V_lv)
plot(times, Vrest)
ylabel('ml')
legend('V lv', 'Vrest')

subplot(5,1,3); 
hold on 
plot(times, P_lv_in_mmHg)
plot(times, P_lv_out_mmHg)
if ~prescribe_Q_ao
    P_ao_distal_mmHg = P_ao_distal / MMHG_TO_CGS; 
    plot(times, P_ao_distal_mmHg); 
end 

ylabel('P mmHg')
legend('P lv in mmHg', 'P lv out mmHg', 'P ao distal mmHg')

subplot(5,1,4); 
hold on 
plot(times, act)
legend('activation')

subplot(5,1,5); 
hold on 
plot(times, Elas)
legend('Elas')




figure; 
hold on 
plot(V_lv, P_lv_in_mmHg)
xlabel('V ml')
ylabel('P LV mmHg')







