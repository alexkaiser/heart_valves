% Script to build valve 


% Size parameter
% Number of points on free edge of each leaflet 
% 
N = 16; 

% Show some output 
plots = true; 

% Initialize structures 
% Many parameters are in this script 

radial = true; 
bead_slip = true; 
attached = false; 

if radial
    
    if bead_slip 
        valve = initialize_valve_data_structures_radial_bead_slip(N, attached); 
    else 
        valve = initialize_valve_data_structures_radial(N); 
    end 
        
else 
    
    if bead_slip || attached 
        error('diagonal fibers not implemented for closed bead slip or attached'); 
    end 
    
    valve = initialize_valve_data_structures(N); 
end 
    
    
if plots 
    fig = surf_plot(valve.posterior); 
    title('Reference configuration of posterior surface'); 
    fig = surf_plot(valve.anterior, fig); 
    title('Reference configuration of anterior surface'); 
end


valve_plot(valve)
if radial
    title('Refernece configuration radial fibers')
else
    title('Refernece configuration diagonal fibers')
end 

iteration_movie_anterior = false; 
if iteration_movie_anterior 
    valve.anterior.iteration_movie = iteration_movie_anterior; 
    valve.anterior.frame = 0; 
    valve.anterior.base_name = sprintf('%s_frames', valve.base_name); 
    valve.anterior.springs_written = false; 
end 

% Can use a scalar pressure 
% Or a range for continuation 
p_range = valve.posterior.p_0; 

% Solve an equilibrium problem for the current X configuration 
if radial && bead_slip && attached
    valve = newton_solve_valve_attached(valve, valve.tol_global, valve.max_it); 
else 
    valve = solve_valve(valve, p_range); 
end 

valve_plot(valve)
if radial
    title('Relaxed configuration radial fibers')
else
    title('Relaxed configuration diagonal fibers')
end 


% Save current data 
save(strcat(valve.base_name, '_final_data')); 

% Write to simulation files 
output_to_ibamr_format(valve); 




