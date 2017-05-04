% Script to build valve 


% Size parameter
% Number of points on free edge of each leaflet 
% 
N = 64; 


% Show some output 
plots = false; 

% Initialize structures 
% Many parameters are in this script 

radial       = true; 
bead_slip    = true; 
attached     = false; 
leaflet_only = false; 
optimization = false; 
decreasing_tension = true; 


if radial
    
    if bead_slip 
        valve = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, decreasing_tension); 
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



if radial && bead_slip && attached
    valve = newton_solve_valve_attached(valve, valve.tol_global, valve.max_it); 
else 
%     valve = solve_valve(valve, p_range, repulsive_coeff_range); 

    [valve valve_with_reference pass_all] = solve_valve(valve); 
end 

fig = figure; 
fig = valve_plot(valve, fig); 

title('Pressurized configuration fibers'); 
saveas(fig, strcat(valve.base_name, '_pressurized'), 'fig'); 

if ~isempty(valve_with_reference)
    fig = figure; 
    fig = valve_plot(valve_with_reference, fig); 
    title('Relaxed configuration radial fibers, reference config based constitutive law'); 
    saveas(fig, strcat(valve.base_name, '_relaxed'), 'fig'); 
end 

if pass_all 
    fprintf('Final solve passed.\n'); 
else 
    fprintf('Final solve failed.\n'); 
end 


% Save current data 
save(strcat(valve.base_name, '_final_data')); 

% Write to simulation files 
if ~isempty(valve_with_reference)
    output_to_ibamr_format(valve_with_reference); 
end 



