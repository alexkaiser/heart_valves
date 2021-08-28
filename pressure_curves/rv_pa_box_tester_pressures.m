
% quadrature spacing 
debug = false; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 
plots = false; 

MMHG_TO_CGS = 1333.22368; 

R_proximal_experimental = 4.89; 
R_distal_experimental   = 46.15; 
C_experimental          = 0.0154; 
R_total_experimental    = R_proximal_experimental + R_distal_experimental; 
tau_experimental        = C_experimental * R_distal_experimental; 

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
bump_radius_rv = .05; 
bump_radius_pa = .05; 
n_fourier_coeffs = 600; 

% for outputting resistance/flow adjusted series 
r_dynespercm2 = 250; 
r_mmHg = r_dynespercm2 / MMHG_TO_CGS; 
r_suffix = "_r_250"; 

points_one_cycle_right_ventricle = [times, rv_pressure]; 
points_one_cycle_pa              = [times, pa_pressure]; 

suffix_right_ventricle = "_right_ventricle"; 
file_name = strcat(base_name, suffix_right_ventricle, '.txt'); 
[a_0_right_ventricle a_n_right_ventricle b_n_right_ventricle Series_right_ventricle times_rv linear_interp_vals_rv] = series_and_smooth(points_one_cycle_right_ventricle, dt, bump_radius_rv, n_fourier_coeffs, plots); 
output_series_coeffs_to_txt(a_0_right_ventricle, a_n_right_ventricle, b_n_right_ventricle, n_fourier_coeffs, cycle_length, file_name); 

suffix_pa = "_pa"; 
file_name = strcat(base_name, suffix_pa, '.txt'); 
[a_0_pa a_n_pa b_n_pa Series_pa times_pa linear_interp_vals_pa] = series_and_smooth(points_one_cycle_pa, dt, bump_radius_pa, n_fourier_coeffs, plots); 
output_series_coeffs_to_txt(a_0_pa, a_n_pa, b_n_pa, n_fourier_coeffs, cycle_length, file_name); 


flows_rv_exp = [-4.679 3.558 34.035 101.834 180.228 226.853 230.468 192.169 123.601 48.604 4.302 -2.023 -0.795 0.197 0.101 -0.209 -1.054 -2.747 -2.837 -3.134]'; 
flows_lpa_exp = [-1.623108 1.62083 13.35831 39.69248 73.04993 94.478 97.37274 82.16487 55.11997 23.08761 1.514264 -2.709484 -2.605704 -3.754506 -4.046469 -4.016213 -4.166175 -4.095049 -3.645273 -2.907113]; 
flows_rpa_exp = [5.198252 6.404848 18.99174 48.68156 84.09315 106.0516 106.0705 87.31557 55.59498 20.70115 1.003064 0.279315 4.22013 5.346351 6.04164 6.95495 7.087463 6.348237 5.817579 5.577208]; 

% puts final elements at beginning 
permute_flows = false; 
if permute_flows
    n_to_permute = 1; 
    flows_rv_exp  = [flows_rv_exp(end-n_to_permute+1:end); flows_rv_exp(1:end-n_to_permute)]; 
    flows_lpa_exp = [flows_lpa_exp(end-n_to_permute+1:end); flows_lpa_exp(1:end-n_to_permute)]; 
    flows_pa_exp  = [flows_rpa_exp(end-n_to_permute+1:end); flows_rpa_exp(1:end-n_to_permute)]; 
end

% adds last flows to time zero 
periodic_wrap_flows = true; 
if permute_flows
    n_to_permute = 1; 
    flows_rv_exp  = [flows_rv_exp(end); flows_rv_exp]; 
    flows_lpa_exp = [flows_lpa_exp(end); flows_lpa_exp]; 
    flows_pa_exp  = [flows_rpa_exp(end); flows_rpa_exp]; 
end

% times_flows = linspace(0, cycle_length, length(flows)); 
times_flows = linspace(cycle_length/length(flows_rv_exp), cycle_length, length(flows_rv_exp)); 

% times are experimental pressure times 
flows_spline_rv  = interp1(times_flows, flows_rv_exp, times, 'spline'); 
flows_spline_rpa = interp1(times_flows, flows_rpa_exp, times, 'spline'); 
flows_spline_lpa = interp1(times_flows, flows_lpa_exp, times, 'spline'); 

rpa_pressure_adjusted_rq = pa_pressure - r_mmHg * flows_spline_rpa; 
lpa_pressure_adjusted_rq = pa_pressure - r_mmHg * flows_spline_lpa; 

points_one_cycle_rpa              = [times, rpa_pressure_adjusted_rq]; 
points_one_cycle_lpa              = [times, lpa_pressure_adjusted_rq]; 

suffix_rpa = "_rpa"; 
file_name = strcat(base_name, suffix_rpa, r_suffix, '.txt'); 
[a_0_rpa a_n_rpa b_n_rpa Series_rpa times_rpa linear_interp_vals_rpa] = series_and_smooth(points_one_cycle_rpa, dt, bump_radius_pa, n_fourier_coeffs, plots); 
output_series_coeffs_to_txt(a_0_rpa, a_n_rpa, b_n_rpa, n_fourier_coeffs, cycle_length, file_name); 

suffix_lpa = "_lpa"; 
file_name = strcat(base_name, suffix_lpa, r_suffix, '.txt'); 
[a_0_lpa a_n_lpa b_n_lpa Series_lpa times_lpa linear_interp_vals_lpa] = series_and_smooth(points_one_cycle_lpa, dt, bump_radius_pa, n_fourier_coeffs, plots); 
output_series_coeffs_to_txt(a_0_lpa, a_n_lpa, b_n_lpa, n_fourier_coeffs, cycle_length, file_name); 



flow_plot = true; 
if flow_plot
    figure; 
    plot(times_flows, flows_rv_exp); 
    hold on; 
    plot(times_flows, flows_rpa_exp); 
    plot(times_flows, flows_lpa_exp); 
    legend('q rv', 'q rpa', 'q lpa')
    title('experimental flows')
end 



basic_series_plots = true; 
if basic_series_plots
    t = 0:dt:cycle_length; 
    vals_right_ventricle_series = Series_right_ventricle(t); 
    fig = figure; 
    hold on 
    plot(t, vals_right_ventricle_series, 'k');

    t = 0:dt:cycle_length; 
    vals_pa_series = Series_pa(t); 
    plot(t, vals_pa_series, 'b'); 
    
    t = 0:dt:cycle_length; 
    vals_rpa_series = Series_rpa(t); 
    plot(t, vals_rpa_series, 'r'); 
    
    t = 0:dt:cycle_length; 
    vals_lpa_series = Series_lpa(t); 
    plot(t, vals_lpa_series, 'g'); 
    
    legend('RV', 'PA', 'RPA adjusted', 'LPA adjusted')
    
    title('RV PA pressure')
    xlabel('t')
    ylabel('p (mmHg)')
    set(fig, 'Position', [100, 100, 1000, 500])
    set(fig,'PaperPositionMode','auto')
    printfig(fig, 'rv_pa_pressure')

%     fig = figure; 
%     plot(t, vals_right_ventricle_series - vals_pa_series, 'k')
%     title('Pressure differences (mmHg)')
end 

interp_flow_plot = false; 
if interp_flow_plot 
    % flows_spline = interp1(times_flows, flows, t, 'spline'); 
    flows_rv_exp_spline = interp1(times_flows, flows_rv_exp, t); 

    figure; 
    plot(times_flows, flows_rv_exp)
    hold on 
    plot(t, flows_rv_exp_spline) 
    title('Flow, experimental trace')
    xlabel('t (s)')
    ylabel('flow (ml/s)')
    legend('pointwise', 'spline')
end 


rcr_estimates = false; 
if rcr_estimates

    dt_flows = cycle_length / length(flows_rv_exp)
    total_flow = sum(flows_rv_exp) * dt_flows % ml 

    Q_mean = total_flow / cycle_length 

    Q_mean_L_per_min = Q_mean * 60/1000 

    ratio_prox_to_distal_resistors = 0.05; % R_proximal_experimental / R_distal_experimental

    decay_time = .4 

    P_mean = mean(vals_pa_series) * MMHG_TO_CGS; 
    
    P_max = 42 * MMHG_TO_CGS; % interpolating the decay by eye to the middle of the oscillation 
                              % in the experimental trace of PA pressure 
    P_min = 31.25 * MMHG_TO_CGS; 

    tol = 1e-12; 

    Q_mean_each = Q_mean / 2; 

    R_coeff = 1.0; 

    % total resistance is determined by mean pressure and mean flow 
    R_total = R_coeff * P_mean / Q_mean_each;  

    % ratio of resistors is constant 
    % resistors sum to total resistance 
    R_distal = R_total / (1.0 + ratio_prox_to_distal_resistors); 
    R_proximal = R_total - R_distal; 

    if abs(R_distal + R_proximal - R_total) > tol 
        error('resistors not adding up correctly')
    end

    C_prefactor = 0.5; 

    % timescale for pressure decrease during aortic valve closure 
    C = C_prefactor * -decay_time / (R_distal * log(P_min/P_max)); 

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
end 

more_plots = false; 
if more_plots

    fig = figure; 

    times_two_cycles = [t,t + cycle_length]; 
    q_rv_exp = [flows_rv_exp_spline, flows_rv_exp_spline]; 

    times = 0:dt:(2*cycle_length);
    p_rv_exp  = Series_right_ventricle(times); 
    p_pa_exp  = Series_pa(times); 

    title('experimental pressures and flows')
    subplot(2,1,1)
    plot(times, p_rv_exp , 'k')
    hold on
    plot(times, p_pa_exp , ':k')

    % interpolant values on top 
    plot(times_rv, linear_interp_vals_rv)
    plot(times_pa, linear_interp_vals_pa)

    legend('P_{RV}', 'P_{PA}', 'rv exp', 'pa exp', 'Location','NorthEastOutside');
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

    % times_exp = times; 

    % this is saving from the series for pressures 
    % save 'bc_variables_experimental.mat' times_two_cycles q_rv_exp times_exp p_rv_exp p_pa_exp

end 






