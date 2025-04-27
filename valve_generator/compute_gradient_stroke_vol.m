function [delta_p_mean, delta_p_max, q_full_cycle, q_systole, q_diff_full_systole, dt, delta_p_mean_alt,delta_p_max_alt] = ...
    compute_gradient_stroke_vol(p_upstream, p_downstream, times_pressure, q_in, times_q, cycle_duration, cycle_number, t_min_pressure, p_downstream_alt)

MMHG_TO_CGS = 1333.22368;

% assumes these are in the current directory 
% integral_metric_data_contour_21; 
% p_lvot = P / (A * MMHG_TO_CGS); 
% 
% integral_metric_data_contour_12; 
% p_st_junction = P / (A * MMHG_TO_CGS); 
% 
% p_diff = p_lvot - p_st_junction; 
% 
% times_pressure = t; 

cycle_start = (cycle_number - 1) * cycle_duration; 
cycle_end   = cycle_number       * cycle_duration; 

t_idx_min = find(times_pressure > t_min_pressure, 1); 

% cycle_start = (cycle_number - 1) * cycle_duration; 
% cycle_end   = cycle_number       * cycle_duration; 

p_diff = p_upstream - p_downstream;


idx_pressure_cross_pos = find((p_diff > 0) & (times_pressure >= times_pressure(t_idx_min)), 1); 
idx_pressure_cross_neg = find((p_diff <= 0) & (times_pressure >= times_pressure(idx_pressure_cross_pos)), 1); 

% 'pressure crossings'
delta_p_max = max(p_diff(idx_pressure_cross_pos:idx_pressure_cross_neg));

% integral on pressure cross times 
delta_p_mean = mean(p_diff(idx_pressure_cross_pos:idx_pressure_cross_neg)); 


if exist('p_downstream_alt', 'var')
    
    p_diff_alt = p_upstream - p_downstream_alt;

    idx_pressure_cross_pos_alt = find((p_diff_alt > 0) & (times_pressure >= times_pressure(t_idx_min)), 1); 
    idx_pressure_cross_neg_alt = find((p_diff_alt <= 0) & (times_pressure >= times_pressure(idx_pressure_cross_pos_alt)), 1); 

    % 'pressure crossings'
    delta_p_max_alt = max(p_diff_alt(idx_pressure_cross_pos_alt:idx_pressure_cross_neg_alt));

    % integral on pressure cross times 
    delta_p_mean_alt = mean(p_diff_alt(idx_pressure_cross_pos_alt:idx_pressure_cross_neg_alt)); 
    
else 
    delta_p_max_alt = nan;
    delta_p_mean_alt = nan;    
end 


min_time_idx_cycle = find(times_q > cycle_start,1);

if cycle_end >= max(times_q)
    max_time_idx_cycle = length(times_q); 
else 
    max_time_idx_cycle = find(times_q >= cycle_end,1); 
end 

q_this_cycle = q_in(min_time_idx_cycle:max_time_idx_cycle);

[max_q max_q_idx] = max(q_this_cycle);


crossing_start_time_found = false;

for idx = max_q_idx:-1:2
    if (q_in(idx-1) <= 0) && (q_in(idx) > 0)
        crossing_start_idx = idx; 
        crossing_start_time_found = true; 
        break;
    end
end 


crossing_end_time_found = false;

for idx = max_q_idx:length(q_this_cycle)
    if (q_in(idx) > 0) && (q_in(idx + 1) < 0)
        crossing_end_idx = idx; 
        crossing_end_time_found = true; 
        break; 
    end
end 

if ~crossing_start_time_found
    error('could not find flow start crossing time'); 
end 

if ~crossing_end_time_found
    error('could not find flow end crossing time'); 
end 


% min_time_idx_cycle = find(times > cycle_start,1);
% 
% if cycle_end >= max(times)
%     max_time_idx_cycle = length(times); 
% else 
%     max_time_idx_cycle = find(times >= cycle_end,1); 
% end 

dt = times_q(2) - times_q(1); 

q_systole = dt * trapz(q_this_cycle(crossing_start_idx:crossing_end_idx)); 

q_full_cycle = dt * trapz(q_in(min_time_idx_cycle:max_time_idx_cycle)); 

q_diff_full_systole = q_full_cycle - q_systole; 


