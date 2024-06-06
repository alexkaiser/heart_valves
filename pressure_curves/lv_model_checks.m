

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
C_max_scaling = 1; 

C_min_ml_over_dynespercm2 = C_min_scaling * C_min_ml_over_kPa * ML_OVER_KPA_TO_ML_OVER_DYNEPERCM2; 
C_max_ml_over_dynespercm2 = C_max_scaling * C_max_ml_over_kPa * ML_OVER_KPA_TO_ML_OVER_DYNEPERCM2; 

Emax = 1/C_max_ml_over_dynespercm2; 
Emin = 1/C_min_ml_over_dynespercm2; 

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

Q_goal_L_per_min = 5.6; 
Q_goal_ml_per_s = Q_goal_L_per_min * 1e3 / 60; 
Q_goal_ml_per_cycle = Q_goal_ml_per_s * cycle_length; 

% quadrature spacing 
debug = true; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 


points_one_cycle_Q_mi = [
0.174135,0.57
0.177188,-15.863
0.180242,-35.697
0.181768,-56.669
0.184821,-73.102
0.187874,-88.4
0.193981,-108.793
0.197034,-87.808
0.200087,-73.627
0.201613,-32.237
0.206193,0.658
0.20772,15.969
0.209246,35.249
0.215352,8.053
0.221459,-10.072
0.229091,6.957
0.24283,-4.911
0.250463,-0.355
0.271835,0.271
0.317631,0.396
0.352742,-0.074
0.392432,13.641
0.410751,42.605
0.426016,94.804
0.430596,131.1
0.439755,169.109
0.453494,204.863
0.485552,171.503
0.519136,131.343
0.551194,94.58
0.583251,66.889
0.619889,49.414
0.673318,34.254
0.735907,27.055
0.786283,23.225
0.812234,30.099
0.819867,45.994
0.8275,71.527
0.833606,101.024
0.835132,128.808
0.845818,155.483
0.850398,173.637
0.86261,150.993
0.86719,120.959
0.868717,94.884
0.879402,63.732
0.883982,31.43
0.890088,-0.868]; 

% start at zero 
times_Q_mi_adjusted = points_one_cycle_Q_mi(:,1) - points_one_cycle_Q_mi(1,1); 

% normalize times to desired cardiac cycle duration 
times_Q_mi_adjusted = times_Q_mi_adjusted * cycle_length / times_Q_mi_adjusted(end); 

t_start_temp = times_Q_mi_adjusted(end); 
t_final_temp = times_Q_mi_adjusted(end); 

Q_mi_adjusted = points_one_cycle_Q_mi(:,2);  

% permute 
% start_idx = find(times_Q_mi_adjusted > .4, 1); 
% 
% first_portion_times_no_adjust = times_Q_mi_adjusted(start_idx:end)
% second_portion_times_no_adjust = times_Q_mi_adjusted(1:start_idx-1)
% 
% first_portion_times = times_Q_mi_adjusted(start_idx:end) % - times_Q_mi_adjusted(start_idx) 
% second_portion_times = times_Q_mi_adjusted(1:start_idx-1) % - times_Q_mi_adjusted(start_idx) 
% 
% % first_portion_times = times_Q_mi_adjusted(start_idx:end) - times_Q_mi_adjusted(start_idx) 
% % second_portion_times = times_Q_mi_adjusted(1:start_idx-1) - times_Q_mi_adjusted(start_idx) + cycle_length 
% times_Q_mi_adjusted = [first_portion_times; second_portion_times]

% Q_mi_adjusted = [Q_mi_adjusted(start_idx:end); Q_mi_adjusted(1:start_idx-1)]


points_one_cycle_Q_mi_adjusted = [times_Q_mi_adjusted, Q_mi_adjusted]; 


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

% scale mitral to goal 

Series_Q_mi_scaled = @(t) (Q_goal_ml_per_cycle / Q_mi_total) * Series_Q_mi(t); 


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

t_offset_Q_mi = -0.25; 

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

prescribe_Q_ao = false; 
if prescribe_Q_ao
    Q_ao = Series_Q_Ao_scaled(times); 
else 
    Q_ao = zeros(n_times, 1); 
    P_ao_distal = Series_P_Ao_dynescm2(times); 
end 

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

    act(n+1) = 0;
    if (t_contract <= t_twitch)
        act(n+1) = (-0.5 * cos(2 * pi * t_contract / t_twitch) + 0.5)^cos_power;
    end 

    Vrest(n+1) = (1.0 - act(n+1)) * (Vrd - Vrs) + Vrs;
    Elas(n+1) = (Emax - Emin) * act(n+1) + Emin;

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
legend('P lv in mmHg', 'P lv out mmHg')

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







