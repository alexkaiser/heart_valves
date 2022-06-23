
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
bump_radius_rv = .01; 
bump_radius_pa = .01;
n_fourier_coeffs = 600; 

% for outputting resistance/flow adjusted series 
r_dynespercm2 = 250; 
r_mmHg = r_dynespercm2 / MMHG_TO_CGS; 
r_suffix = "_r_250"; 

% IB flow lost 
radius_valve = 1; 
dx_fluid = 0.045 * 2; 
frac_flow_expected = (radius_valve - dx_fluid)^4; 
coeff_pressure_diff_adjust = radius_valve^4 / frac_flow_expected - 1; 

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

% points_one_cycle_difference = [times, rv_pressure - pa_pressure]; 
% [a_0_difference a_n_difference b_n_difference Series_difference times_difference linear_interp_vals_difference] = series_and_smooth(points_one_cycle_difference, dt, bump_radius_pa, n_fourier_coeffs, plots); 

% smoothed differences 
pressure_diff_positive = max(rv_pressure - pa_pressure, 0); 
points_one_cycle_difference_positive = [times, pressure_diff_positive]; 
[a_0_difference_positive a_n_difference_positive b_n_difference_positive Series_difference_positive times_difference_positive linear_interp_vals_differenc_positivee] = series_and_smooth(points_one_cycle_difference_positive, dt, bump_radius_pa, n_fourier_coeffs, plots); 

pressure_diff = rv_pressure - pa_pressure; 
points_one_cycle_difference = [times, pressure_diff]; 
[a_0_difference a_n_difference b_n_difference Series_difference times_difference linear_interp_vals_differenc] = series_and_smooth(points_one_cycle_difference, dt, bump_radius_pa, n_fourier_coeffs, plots); 


% rv adjusted 
pressure_rv_plus_diff_positive = rv_pressure + coeff_pressure_diff_adjust * pressure_diff_positive; 
points_one_cycle_right_ventricle_plus_diff_positive = [times, pressure_rv_plus_diff_positive]; 
suffix_right_ventricle__plus_diff_positive = "_right_ventricle_plus_diff_positive"; 
file_name = strcat(base_name, suffix_right_ventricle__plus_diff_positive, '.txt'); 
[a_0_right_ventricle_plus_diff_positive a_n_right_ventricle_plus_diff_positive  b_n_right_ventricle_plus_diff_positive  ...
    Series_right_ventricle_plus_diff_positive  times_rv_plus_diff_positive  linear_interp_vals_rv_plus_diff_positive] = ...
    series_and_smooth(points_one_cycle_right_ventricle_plus_diff_positive, dt, bump_radius_rv, n_fourier_coeffs, plots); 
output_series_coeffs_to_txt(a_0_right_ventricle_plus_diff_positive, a_n_right_ventricle_plus_diff_positive, b_n_right_ventricle_plus_diff_positive, n_fourier_coeffs, cycle_length, file_name); 





flows_rv_exp = [-4.679 3.558 34.035 101.834 180.228 226.853 230.468 192.169 123.601 48.604 4.302 -2.023 -0.795 0.197 0.101 -0.209 -1.054 -2.747 -2.837 -3.134]'; 
flows_lpa_exp = [-1.623108 1.62083 13.35831 39.69248 73.04993 94.478 97.37274 82.16487 55.11997 23.08761 1.514264 -2.709484 -2.605704 -3.754506 -4.046469 -4.016213 -4.166175 -4.095049 -3.645273 -2.907113]'; 
flows_rpa_exp = [5.198252 6.404848 18.99174 48.68156 84.09315 106.0516 106.0705 87.31557 55.59498 20.70115 1.003064 0.279315 4.22013 5.346351 6.04164 6.95495 7.087463 6.348237 5.817579 5.577208]'; 




% puts final elements at beginning 
permute_flows = false; 
if permute_flows
    n_to_permute = 1; 
    flows_rv_exp  = [flows_rv_exp(end-n_to_permute+1:end); flows_rv_exp(1:end-n_to_permute)]; 
    flows_lpa_exp = [flows_lpa_exp(end-n_to_permute+1:end); flows_lpa_exp(1:end-n_to_permute)]; 
    flows_pa_exp  = [flows_rpa_exp(end-n_to_permute+1:end); flows_rpa_exp(1:end-n_to_permute)]; 
end

% % adds last flows to time zero 
% periodic_wrap_flows = true; 
% if periodic_wrap_flows
%     flows_rv_exp  = [flows_rv_exp(end); flows_rv_exp]; 
%     flows_lpa_exp = [flows_lpa_exp(end); flows_lpa_exp]; 
%     flows_pa_exp  = [flows_rpa_exp(end); flows_rpa_exp]; 
% end

% times_flows = linspace(0, cycle_length, length(flows)); 
times_flows = linspace(cycle_length/length(flows_rv_exp), cycle_length, length(flows_rv_exp)); 



    
% times are experimental pressure times 
flows_spline_rv  = interp1(times_flows, flows_rv_exp, times, 'spline'); 
flows_spline_rpa = interp1(times_flows, flows_rpa_exp, times, 'spline'); 
flows_spline_lpa = interp1(times_flows, flows_lpa_exp, times, 'spline'); 


resistance_adjusted_pressures = false; 

if resistance_adjusted_pressures

    individual_flow_adjust = false; 
    if individual_flow_adjust
        rpa_pressure_adjusted_rq = pa_pressure - r_mmHg * flows_spline_rpa; 
        lpa_pressure_adjusted_rq = pa_pressure - r_mmHg * flows_spline_lpa; 
    else 
        rpa_pressure_adjusted_rq = pa_pressure - r_mmHg * flows_spline_rv/2; 
        lpa_pressure_adjusted_rq = pa_pressure - r_mmHg * flows_spline_rv/2; 
    end 

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
end 


% poiseuille flow estimates with paraview hack radii 
% do not use these 
poiseuille_flow_est = false; 
% poiseuille_adjusted_pressures = false; 
% if poiseuille_flow_est
% 
%     mu = 0.039; 
%     
%     q_max_rpa = max(flows_rpa_exp); 
%     q_max_lpa = max(flows_lpa_exp); 
%     
%     len_lpa = 7.71; 
%     rad_lpa = 1.35/2; 
%     
%     len_rpa = 3.08; 
%     rad_rpa = 1.47/2; 
%     
%     deltap_at_qmax_dynes_rpa = 8 * mu * len_rpa * q_max_rpa / (pi * rad_rpa^4); 
%     deltap_at_qmax_rpa = deltap_at_qmax_dynes_rpa / MMHG_TO_CGS 
%     
%     deltap_at_qmax_dynes_lpa = 8 * mu * len_lpa * q_max_lpa / (pi * rad_lpa^4);     
%     deltap_at_qmax_lpa = deltap_at_qmax_dynes_lpa / MMHG_TO_CGS
%     
%     rad_lpa_effective_ib_adjust = rad_lpa - dx_fluid; 
%     
%     deltap_at_qmax_dynes_ib_adjust_lpa = 8 * mu * len_lpa * q_max_lpa / (pi * rad_lpa_effective_ib_adjust^4);     
%     deltap_at_qmax_ib_adjust_lpa = deltap_at_qmax_dynes_ib_adjust_lpa / MMHG_TO_CGS
% 
%        
%     if poiseuille_adjusted_pressures
%         
%         rpa_pressure_adjusted_rq = pa_pressure - deltap_at_qmax_rpa; 
%         lpa_pressure_adjusted_rq = pa_pressure - deltap_at_qmax_lpa; 
% 
%         points_one_cycle_rpa              = [times, rpa_pressure_adjusted_rq]; 
%         points_one_cycle_lpa              = [times, lpa_pressure_adjusted_rq]; 
%         
%         suffix_rpa = "_rpa"; 
%         file_name = strcat(base_name, suffix_rpa, '.txt'); 
%         [a_0_rpa a_n_rpa b_n_rpa Series_rpa times_rpa linear_interp_vals_rpa] = series_and_smooth(points_one_cycle_rpa, dt, bump_radius_pa, n_fourier_coeffs, plots); 
%         output_series_coeffs_to_txt(a_0_rpa, a_n_rpa, b_n_rpa, n_fourier_coeffs, cycle_length, file_name); 
% 
%         suffix_lpa = "_lpa"; 
%         file_name = strcat(base_name, suffix_lpa, '.txt'); 
%         [a_0_lpa a_n_lpa b_n_lpa Series_lpa times_lpa linear_interp_vals_lpa] = series_and_smooth(points_one_cycle_lpa, dt, bump_radius_pa, n_fourier_coeffs, plots); 
%         output_series_coeffs_to_txt(a_0_lpa, a_n_lpa, b_n_lpa, n_fourier_coeffs, cycle_length, file_name); 
%     end 
%     
% end 


% IB estimates at changes of linear resistance 
ib_linear_resistance_estimates = true; 
poiseuille_adjusted_pressures = true;
rv_adjusted = false; 
if ib_linear_resistance_estimates  
    
    mu = 0.039; 
    
    % compute radius from known resistance and length 
    radius_poiseuille = @(L, resistance) (8*mu*L/(pi*resistance))^(1/4); 
    
    delta_pressure_poiseuille = @(L,q,radius) 8 * mu * L * q / (pi * radius^4);
    
    resistance_poiseuille = @(L,radius) 8 * mu * L / (pi * radius^4);
    
    
    % basic resistances from SV 
%     L_rv = 14.759353372237321; 
%     resistance_rv = 0.09405289588937868; 

    % redone RV resistances with valve area cropped out 
    L_rv = 10.158089589432105; 
    resistance_rv = 0.05603465074128439; 
    r_rv = radius_poiseuille(L_rv, resistance_rv) 
    
    L_valve = 1.5; % total scaffold height 
    r_valve = 1.0; 
    resistance_valve = resistance_poiseuille(L_valve,r_valve); 
    
    L_rpa = 4.519774393313929; 
    resistance_rpa = 1.6645415607251934; 
    r_rpa = radius_poiseuille(L_rpa, resistance_rpa) 
    
    L_lpa = 8.139681927111639; 
    resistance_lpa = 3.877335669649555; 
    r_lpa = radius_poiseuille(L_lpa, resistance_lpa) 

    % pressure diffs at nominal flow rates 
    q_mean_nominal_l_per_min = 3.5; 
    q_mean_nominal = q_mean_nominal_l_per_min * 1000 / 60; 
    
    deltap_at_q_nominal_dynes_rpa = delta_pressure_poiseuille(L_rpa, q_mean_nominal, r_rpa);
    deltap_at_q_nominal_rpa = deltap_at_q_nominal_dynes_rpa / MMHG_TO_CGS 
    
    deltap_at_q_nominal_dynes_lpa = delta_pressure_poiseuille(L_lpa, q_mean_nominal, r_lpa);
    deltap_at_q_nominal_lpa = deltap_at_q_nominal_dynes_lpa / MMHG_TO_CGS
    

    % compute systole/diastole times 
    systole_start = fzero(Series_difference, 0.1);
    systole_end = fzero(Series_difference, 0.4);
    systole_duration = systole_end - systole_start
    
    systole_time_fraction = systole_duration / cycle_length 
    
    % expected pressure differences with Poiseuille flow in systole only 
    deltap_at_q_nominal_systole_rpa = deltap_at_q_nominal_rpa / systole_time_fraction
    deltap_at_q_nominal_systole_lpa = deltap_at_q_nominal_lpa / systole_time_fraction
    
    % check answers 
%     resistance_rpa_formula = resistance_poiseuille(L_rpa, r_rpa)
%     resistance_lpa_formula = resistance_poiseuille(L_lpa, r_lpa)
    
    
    % increase estimates by IB resistance ratios 
    resistance_IB_adjust_rv    = resistance_poiseuille(L_rv, r_rv - dx_fluid)
    resistance_IB_adjust_valve = resistance_poiseuille(L_valve, r_valve - dx_fluid)
    resistance_IB_adjust_rpa   = resistance_poiseuille(L_rpa, r_rpa - dx_fluid)
    resistance_IB_adjust_lpa   = resistance_poiseuille(L_lpa, r_lpa - dx_fluid)
    
    % pressure differences with resistance increases 
    delta_p_q_nominal_systole_IB_adjust_rpa = (resistance_IB_adjust_rpa/resistance_rpa) * deltap_at_q_nominal_systole_rpa 
    delta_p_q_nominal_systole_IB_adjust_lpa = (resistance_IB_adjust_lpa/resistance_lpa) * deltap_at_q_nominal_systole_lpa 
    
    reistance_ratio_rv_valve = (resistance_IB_adjust_rv + resistance_IB_adjust_valve) / (resistance_rv + resistance_valve)
    
    
    if poiseuille_adjusted_pressures
        
        mean_pressure_diff = mean(pressure_diff_positive(pressure_diff_positive > 0)); 
        
        pressure_diff_positive_normalized = pressure_diff_positive / mean_pressure_diff; 
        
        rpa_pressure_adjusted_q_nomincal_systole = pa_pressure - deltap_at_q_nominal_systole_rpa * pressure_diff_positive_normalized; 
        lpa_pressure_adjusted_r_nomincal_systole = pa_pressure - deltap_at_q_nominal_systole_lpa * pressure_diff_positive_normalized; 

        points_one_cycle_rpa              = [times, rpa_pressure_adjusted_q_nomincal_systole]; 
        points_one_cycle_lpa              = [times, lpa_pressure_adjusted_r_nomincal_systole]; 

        suffix_rpa = "_rpa"; 
        file_name = strcat(base_name, suffix_rpa, '.txt'); 
        [a_0_rpa a_n_rpa b_n_rpa Series_rpa times_rpa linear_interp_vals_rpa] = series_and_smooth(points_one_cycle_rpa, dt, bump_radius_pa, n_fourier_coeffs, plots); 
        output_series_coeffs_to_txt(a_0_rpa, a_n_rpa, b_n_rpa, n_fourier_coeffs, cycle_length, file_name); 

        suffix_lpa = "_lpa"; 
        file_name = strcat(base_name, suffix_lpa, '.txt'); 
        [a_0_lpa a_n_lpa b_n_lpa Series_lpa times_lpa linear_interp_vals_lpa] = series_and_smooth(points_one_cycle_lpa, dt, bump_radius_pa, n_fourier_coeffs, plots); 
        output_series_coeffs_to_txt(a_0_lpa, a_n_lpa, b_n_lpa, n_fourier_coeffs, cycle_length, file_name);         
    end 
    
    
    if rv_adjusted
        coeff_pressure_diff_adjust = reistance_ratio_rv_valve - 1; 
        pressure_rv_plus_diff_positive = rv_pressure + coeff_pressure_diff_adjust * pressure_diff_positive; 
        points_one_cycle_right_ventricle_plus_diff_positive = [times, pressure_rv_plus_diff_positive]; 
        suffix_right_ventricle__plus_diff_positive = "_right_ventricle_rv_valve_ib_adjust"; 
        file_name = strcat(base_name, suffix_right_ventricle__plus_diff_positive, '.txt'); 
        [a_0_right_ventricle_plus_diff_positive a_n_right_ventricle_plus_diff_positive  b_n_right_ventricle_plus_diff_positive  ...
            Series_right_ventricle_plus_diff_positive  times_rv_plus_diff_positive  linear_interp_vals_rv_plus_diff_positive] = ...
            series_and_smooth(points_one_cycle_right_ventricle_plus_diff_positive, dt, bump_radius_rv, n_fourier_coeffs, plots); 
        output_series_coeffs_to_txt(a_0_right_ventricle_plus_diff_positive, a_n_right_ventricle_plus_diff_positive, b_n_right_ventricle_plus_diff_positive, n_fourier_coeffs, cycle_length, file_name); 
    end 
    
    
    reduced_radius_model_estimates = false; 
    if reduced_radius_model_estimates
        % resistances from SV on reduced model radius of 0.09 cm 
        L_rv_reduced_pt09cm = 12.795778735700885; 
        resistance_rv_reduced_pt09cm = 0.13410025248596413; 
        r_rv_reduced_pt09cm = radius_poiseuille(L_rv_reduced_pt09cm , resistance_rv_reduced_pt09cm ) 

        L_rpa_reduced_pt09cm = 4.526318740007514; 
        resistance_rpa_reduced_pt09cm = 2.9466135041324772; 
        r_rpa_reduced_pt09cm  = radius_poiseuille(L_rpa_reduced_pt09cm , resistance_rpa_reduced_pt09cm ) 

        L_lpa_reduced_pt09cm = 8.130325343022083; 
        resistance_lpa_reduced_pt09cm = 7.1320202349241235; 
        r_lpa_reduced_pt09cm = radius_poiseuille(L_lpa_reduced_pt09cm , resistance_lpa_reduced_pt09cm ) 

        r_rv_diff = r_rv - r_rpa_reduced_pt09cm 
        r_lpa_diff = r_lpa - r_lpa_reduced_pt09cm 
        r_rpa_diff = r_rpa - r_rpa_reduced_pt09cm

        resistance_ratio_rv = resistance_rv_reduced_pt09cm / resistance_rv 
        resistance_ratio_lpa = resistance_lpa_reduced_pt09cm / resistance_lpa
        resistance_ratio_rpa = resistance_rpa_reduced_pt09cm / resistance_rpa 

        q_max_rv  = max(flows_rv_exp); 
        q_max_lpa = max(flows_lpa_exp); 
        q_max_rpa = max(flows_rpa_exp); 

        deltap_rv_sv_resistance = resistance_rv * q_max_rv;
        deltap_lpa_sv_resistance = resistance_lpa * q_max_lpa;
        deltap_rpa_sv_resistance = resistance_rpa * q_max_rpa;
        deltap_rv_sv_resistance_mmHg = deltap_rv_sv_resistance / MMHG_TO_CGS
        deltap_lpa_sv_resistance_mmHg = deltap_lpa_sv_resistance / MMHG_TO_CGS
        deltap_rpa_sv_resistance_mmHg = deltap_rpa_sv_resistance / MMHG_TO_CGS

        deltap_IB_rv_sv_resistance = resistance_ratio_lpa * deltap_rv_sv_resistance;
        deltap_IB_lpa_sv_resistance = resistance_ratio_lpa * deltap_lpa_sv_resistance;
        deltap_IB_rpa_sv_resistance = resistance_ratio_rpa * deltap_rpa_sv_resistance;
        deltap_IB_rv_sv_resistance_mmHg = deltap_IB_rv_sv_resistance / MMHG_TO_CGS 
        deltap_IB_lpa_sv_resistance_mmHg = deltap_IB_lpa_sv_resistance / MMHG_TO_CGS 
        deltap_IB_rpa_sv_resistance_mmHg = deltap_IB_rpa_sv_resistance / MMHG_TO_CGS 
    end 
end 


if resistance_adjusted_pressures && poiseuille_adjusted_pressures
    error('resistance_adjusted_pressures && poiseuille_adjusted_pressures cannot both be true')
end 

flow_plot = true; 
if flow_plot
    fig = figure; 
    plot(times_flows, flows_rv_exp); 
    hold on; 
    plot(times_flows, flows_rpa_exp); 
    plot(times_flows, flows_lpa_exp); 
    
    total_flow = flows_rv_exp - flows_rpa_exp - flows_lpa_exp; 
    plot(times_flows, total_flow)
    
    plot(times_flows, flows_rv_exp/2); 
    
    xlim([0 .835])
    legend('q rv', 'q rpa', 'q lpa', 'total', 'qrv/2')
    title('experimental flows')
    printfig(fig, 'exp_flow_pa_4dmri')
end 

pressure_plot = true; 
if pressure_plot
    fig = figure; 
    times_exp = [table.Time]; 
    p_rv_exp = [table.RightVentriclePressure_Inlet_]; 
    p_pa_exp = [table.MainPulmonaryArteryPressure_Outlets_];
    
    plot(times_exp, p_rv_exp)
    hold on 
    plot(times_exp, p_pa_exp)
    xlim([0 .835])
    legend('P RV', 'P PA')
    title('experimental pressures')
    printfig(fig, 'exp_pressure_4dmri')
    
    fig = figure; 
    plot(times_exp, p_rv_exp - p_pa_exp); 
    axis([0 .4 0 20])    
    title('experimental pressure difference')
    max_fwd_delta_p = max(p_rv_exp - p_pa_exp)
    
    
    
end 

output_experimental_flows = true; 
if output_experimental_flows    
    % output the data directly instead 
    times_exp = [table.Time; table.Time + table.Time(end)]; 
    p_rv_exp = [table.RightVentriclePressure_Inlet_; table.RightVentriclePressure_Inlet_]; 
    p_pa_exp = [table.MainPulmonaryArteryPressure_Outlets_; table.MainPulmonaryArteryPressure_Outlets_]; 

    times_two_cycles = [times_flows, times_flows + cycle_length]'; 
    q_rv_exp = [flows_rv_exp; flows_rv_exp];  
    q_rpa_exp = [flows_rpa_exp; flows_rpa_exp];  
    q_lpa_exp = [flows_lpa_exp; flows_lpa_exp];  
        
    save 'bc_variables_experimental.mat' times_two_cycles q_rv_exp q_rpa_exp q_lpa_exp times_exp p_rv_exp p_pa_exp
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
    
    if resistance_adjusted_pressures || poiseuille_adjusted_pressures
        t = 0:dt:cycle_length; 
        vals_rpa_series = Series_rpa(t); 
        plot(t, vals_rpa_series, 'r'); 

        t = 0:dt:cycle_length; 
        vals_lpa_series = Series_lpa(t); 
        plot(t, vals_lpa_series, 'g'); 
    end 
    
%     t = 0:dt:cycle_length; 
%     vals_rv_plus_diff_series = Series_right_ventricle(t) + coeff_pressure_diff_adjust * Series_difference(t); 
%     plot(t, vals_rv_plus_diff_series); 
%     
    t = 0:dt:cycle_length; 
    vals_rv_plus_diff_positive_series = Series_right_ventricle_plus_diff_positive(t); 
    plot(t, vals_rv_plus_diff_positive_series); 
    
    times_exp = [table.Time]; 
    p_rv_exp = [table.RightVentriclePressure_Inlet_]; 
    p_pa_exp = [table.MainPulmonaryArteryPressure_Outlets_];    
    plot(times_exp, p_rv_exp)
    plot(times_exp, p_pa_exp)
    
    if resistance_adjusted_pressures || poiseuille_adjusted_pressures
        legend('RV', 'PA', 'RPA adjusted', 'LPA adjusted', 'RV plus diff postive', 'P RV EXP', 'P PA EXP')
    else 
        legend('RV', 'PA', 'RV plus diff postive', 'P RV EXP', 'P PA EXP')
    end 
    
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






