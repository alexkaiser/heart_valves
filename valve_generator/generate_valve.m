% Script to build valve 


% Size parameter
% Number of points on free edge of each leaflet 
% 
N = 128; 

% Show some output 
plots = false; 

% Initialize structures 
% Many parameters are in this script 

radial       = true; 
bead_slip    = true; 
attached     = false; 
leaflet_only = false; 
optimization = false; 
repulsive_potential = false; 
decreasing_tension = true; 


if radial
    
    if bead_slip 
        valve = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, repulsive_potential, decreasing_tension); 
    else        
        valve = initialize_valve_data_structures_radial(N); 
    end 
        
else 
    
    if bead_slip || attached || leaflet_only 
        error('diagonal fibers not implemented for closed bead slip or attached'); 
    end 
    
    valve = initialize_valve_data_structures(N); 
end 
    
    
if plots 
    fig = surf_plot(valve.posterior); 
    title('Reference configuration of posterior surface'); 
    fig = surf_plot(valve.anterior, fig); 
    title('Reference configuration of anterior surface'); 
    
    valve_plot(valve)
    if radial
        title('Refernece configuration radial fibers')
    else
        title('Refernece configuration diagonal fibers')
    end 
    
end


iteration_movie_anterior = false; 
if iteration_movie_anterior 
    valve.anterior.iteration_movie = iteration_movie_anterior; 
    valve.anterior.frame = 0; 
    valve.anterior.base_name = sprintf('%s_frames', valve.base_name); 
    valve.anterior.springs_written = false; 
end 


load mitral_tree_128_final_data; 

% Can use a scalar pressure 
% Or a range for continuation 
p_range = valve.anterior.p_0 .* [0:.1:.9, .925:.025:1]; 
% p_range = valve.posterior.p_0; 

linear_open_config  = true; 
p_range_linear      = valve.anterior.p_0; % .* (1:-.1:0); 
strain = 0.16; 
repulsive_coeff_range = []; % [.9:(-0.1):.1]; 


if radial && bead_slip && attached
    valve = newton_solve_valve_attached(valve, valve.tol_global, valve.max_it); 
else 
%     valve = solve_valve(valve, p_range, repulsive_coeff_range); 

    [valve valve_linear pass_all] = solve_valve(valve, p_range, linear_open_config, p_range_linear, strain); 
end 

fig = figure; 
fig = valve_plot(valve, fig); 
if radial
    title('Pressurized configuration radial fibers')
else
    title('Relaxed configuration diagonal fibers')
end 

fig = figure; 
fig = valve_plot(valve_linear, fig); 
title('Relaxed configuration radial fibers, linear constitutive law'); 


if pass_all 
    fprintf('Final solve passed.\n'); 
else 
    fprintf('Final solve failed.\n'); 
end 


% Save current data 
save(strcat(valve.base_name, '_final_data')); 

% Write to simulation files 
output_to_ibamr_format(valve_linear); 




