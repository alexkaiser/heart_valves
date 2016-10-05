% Script to build valve 


% Size parameter
% Number of points on free edge of each leaflet 
% 
N = 32; 

% Show some output 
plots = true; 

% Initialize structures 
% Many parameters are in this script 

radial = true; 
closed_bead_slip = false; 

if radial
    
    if closed_bead_slip 
        valve = initialize_valve_data_structures_radial_bead_slip(N); 
    else 
        valve = initialize_valve_data_structures_radial(N); 
    end 
        
else 
    
    if closed_bead_slip 
        error('diagonal fibers not implemented for closed bead slip'); 
    end 
    
    valve = initialize_valve_data_structures(N); 
end 
    
    
if plots 
    fig = surf_plot(valve.posterior); 
    title('Reference configuration of posterior surface'); 
    fig = surf_plot(valve.anterior); 
    title('Reference configuration of anterior surface'); 
end
    
% Can use a scalar pressure 
% Or a range for continuation 
p_range = valve.posterior.p_0; 

% Solve an equilibrium problem for the current X configuration 
valve = solve_valve(valve, p_range); 

% Save current data 
save(strcat(valve.base_name, '_final_data')); 

% Write to simulation files 
output_to_ibamr_format(valve); 




