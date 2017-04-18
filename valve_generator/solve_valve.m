function [valve valve_with_reference pass_all] = solve_valve(valve, strain)
% 
% Refines valve data structure to equilibrium 
% Applies auto-continuation to pressure and updates both leaflets 
% 

pass_all = true; 


for i=1:length(valve.leaflets)
    
    leaflet = valve.leaflets(i); 
    
    p_initial = 0; 
    p_goal    = leaflet.p_0; 

    [valve.leaflets(i) pass err] = solve_valve_pressure_auto_continuation(leaflet, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

    if pass
        fprintf('Global solve passed, err = %e\n\n', err); 
    else 
        fprintf('Global solve failed, err = %e\n\n', err); 
    end 
    
    fig = figure; 
    surf_plot(valve.leaflets(i), fig); 
    pause(0.01);
    
    pass_all = pass_all && pass; 
    
end 

% constitutive law version 
valve_with_reference = valve; 

% kill off the old leaflet structure, new one has different fields, 
% which makes matlab complain about assigning it to a structure array 
valve_with_reference = rmfield(valve_with_reference, 'leaflets'); 

for i=1:length(valve.leaflets)
        
    valve_with_reference.leaflets(i) = set_rest_lengths_and_constants(valve.leaflets(i), strain); 
    
    leaflet = valve_with_reference.leaflets(i); 
    
    p_initial = leaflet.p_0; 
    p_goal    = 0;

    [valve_with_reference.leaflets(i) pass err] = solve_valve_pressure_auto_continuation(leaflet, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

    if pass
        fprintf('Global solve passed, err = %e\n\n', err); 
    else 
        fprintf('Global solve failed, err = %e\n\n', err); 
    end 
    
    fig = figure; 
    surf_plot(valve.leaflets(i), fig); 
    pause(0.01);
    
    pass_all = pass_all && pass; 
    
end 


