function [delta_p_mean, delta_p_max, q_full_cycle, q_systole, q_diff_full_systole, dt] = compute_gradient_stroke_vol(cycle_duration, cycle_number, t_min_pressure, t_max_pressure, t_min_flow, t_mid_flow)

MMHG_TO_CGS = 1333.22368;

% assumes these are in the current directory 
integral_metric_data_contour_21; 
p_lvot = P / (A * MMHG_TO_CGS); 

integral_metric_data_contour_12; 
p_st_junction = P / (A * MMHG_TO_CGS); 

p_diff = p_lvot - p_st_junction; 

times_pressure = t; 

% mid systole ranges 
t_idx_min = find(times_pressure > t_min_pressure, 1); 
t_idx_max = find(times_pressure > t_max_pressure, 1); 

cycle_start = (cycle_number - 1) * cycle_duration; 
cycle_end   = cycle_number       * cycle_duration; 

t_idx_min_cycle = find(times_pressure > cycle_start,1); 

if cycle_end >= max(times_pressure)
    t_idx_max_cycle = length(times_pressure); 
else 
    t_idx_max_cycle = find(times_pressure > cycle_end,1); 
end 

delta_p_max = max(p_diff(t_idx_min_cycle:t_idx_max_cycle)); 

delta_p_mean = mean(p_diff(t_idx_min:t_idx_max)); 

% flows 
load('bc_data.mat', 'times', 'q_aorta'); 

min_time_idx = find(times > t_min_flow, 1); 
% max_time_idx = find(times > t_max, 1); 

mid_time_idx = find(times > t_mid_flow,1); 

crossing_time_found = false;

for idx = mid_time_idx:(length(times)-1)
    if (q_aorta(idx) > 0) && (q_aorta(idx + 1) < 0)
        crossing_idx = idx; 
        crossing_time_found = true; 
    end 
end 

if ~crossing_time_found
    error('could not find flow crossing time'); 
end 


min_time_idx_cycle = find(times > cycle_start,1);

if cycle_end >= max(times)
    max_time_idx_cycle = length(times); 
else 
    max_time_idx_cycle = find(times >= cycle_end,1); 
end 

dt = times(2) - times(1); 

q_systole = dt * trapz(q_aorta(min_time_idx:crossing_idx)); 

q_full_cycle = dt * trapz(q_aorta(min_time_idx_cycle:max_time_idx_cycle)); 

q_diff_full_systole = q_full_cycle - q_systole; 
