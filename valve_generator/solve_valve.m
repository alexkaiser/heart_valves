function [valve valve_with_reference pass_all] = solve_valve(valve)
% 
% Refines valve data structure to equilibrium 
% Applies auto-continuation to pressure and updates both leaflets 
% 

pass_all = true; 

if length(valve.leaflets) ~= 1
    error('running with single leaflet assumption for now'); 
end 

tol_global            = valve.tol_global; 
max_it                = valve.max_it; 
max_it_continuation   = valve.max_it_continuation; 
max_consecutive_fails = valve.max_consecutive_fails; 
max_total_fails       = valve.max_total_fails; 



for i=1:length(valve.leaflets)
    
    leaflet = valve.leaflets(i); 
    
    p_initial = 0; 
    p_goal    = leaflet.p_0; 

    [valve.leaflets(i) pass err] = solve_valve_pressure_auto_continuation(leaflet, tol_global, max_it, max_it_continuation, p_initial, p_goal, max_consecutive_fails, max_total_fails); 

    if pass
        fprintf('Global solve passed, err = %e\n\n', err); 
    else 
        fprintf('Global solve failed, err = %e\n\n', err); 
    end 
    
%     fig = figure; 
%     surf_plot(valve.leaflets(i), fig); 
%     pause(0.01);
    
    pass_all = pass_all && pass; 
    
end 

if valve.interactive && pass_all 
    fprintf('Solve passed, interactive mode enabled.'); 
    
    fig = figure; 
    valve_plot(valve, fig); 
    title('Valve in interactive mode'); 
    
    while true 
    
        fprintf('Current tension_coeffs struct, which includes all valid variables:\n')
        valve.leaflets(1).tension_coeffs

         try 
            var_name = input('Enter the name of variable to change as a string (no quotes or whitespace, must match exactly):\n', 's'); 

            if ~ischar(var_name)
                fprintf('Must input a valid string for variable name.\n'); 
                continue; 
            end 

            if isempty(var_name)
                var_name = input('No variable name entered. Enter empty string again to leave interactive mode.\n', 's'); 
                if isempty(var_name)
                    break; 
                end 
            end 

            if isfield(valve.leaflets(1).tension_coeffs, var_name)

                tension_coeffs_current = valve.leaflets(1).tension_coeffs; 
                
                if length(tension_coeffs_current.(var_name)) > 1
                    fprintf('Found valid variable %s. You have selected an array variable. Contents of array are:', var_name);
                    tension_coeffs_current.(var_name)
                    idx = input('Enter the index you would like to change, ranges okay:\n');  
                    value = input('Input new value:\n'); 
                    tension_coeffs_current.(var_name)(idx) = value;
                else                    
                    value_old = tension_coeffs_current.(var_name);
                    fprintf('Found valid variable %s with old value %f.\n', var_name, value_old); 
                    value = input('Input new value:\n'); 
                    tension_coeffs_current.(var_name) = value;
                end 
                
                [leaflet_current valve_current] = set_tension_coeffs(valve.leaflets(1), valve, tension_coeffs_current); 

                try
                    [leaflet_current pass err] = newton_solve_valve(leaflet_current, tol_global, max_it, max_consecutive_fails, max_total_fails);  

                    if pass 
                        valve             = valve_current; 
                        valve.leaflets(1) = leaflet_current; 
                        [az el] = view; 
                        clf(fig); 
                        valve_plot(valve, fig);
                        view(az,el); 
                    else 
                        fprintf('New parameters failed. Keeping old tension structure. No pass flag.\n'); 
                    end
                catch 
                    fprintf('New parameters failed. Keeping old tension structure. Catch block, error called in Newton solve.\n'); 
                end 

            else 
                fprintf('Variable not found.\n'); 
            end 
            
         catch 
             fprintf('Error in interactive loop of some kind, restart loop.\n'); 
         end 
    end 
end 



% constitutive law version 
valve_with_reference = valve; 

% kill off the old leaflet structure, new one has different fields, 
% which makes matlab complain about assigning it to a structure array 
valve_with_reference = rmfield(valve_with_reference, 'leaflets'); 

for i=1:length(valve.leaflets)
        
    valve_with_reference.leaflets(i) = set_rest_lengths_and_constants(valve.leaflets(i), valve); 
    
    leaflet = valve_with_reference.leaflets(i); 
    
    p_initial = leaflet.p_0/5; 
    p_goal    = leaflet.p_0/100; 

    [valve_with_reference.leaflets(i) pass err] = solve_valve_pressure_auto_continuation(leaflet, tol_global, max_it, max_it_continuation, p_initial, p_goal, max_consecutive_fails, max_total_fails); 

    if pass
        fprintf('Global solve passed, err = %e\n\n', err); 
    else 
        fprintf('Global solve failed, err = %e\n\n', err); 
    end 
    
%     fig = figure; 
%     surf_plot(valve.leaflets(i), fig); 
%     pause(0.01);
    
    pass_all = pass_all && pass; 
    
end 


