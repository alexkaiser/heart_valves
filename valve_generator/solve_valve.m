function valve = solve_valve(valve, p_range)
% 
% Refines valve data structure to equilibrium 
% Applies auto-continuation to ref_frac and updates both leaflets 
% 

for p_0 = p_range

    % when the pressure changes, just update the pressure and re-run the setup 
    valve.posterior.p_0 = p_0; 
    valve.posterior.p_0 = p_0; 

    [valve.posterior pass err_posterior] = solve_valve_auto_continuation(valve.posterior, valve.tol_global, 'posterior'); 

    if pass 
        fprintf('Global solve passed posterior, err = %f\n\n', err_posterior); 
    else 
        fprintf('Global solve failed posterior, err = %f\n\n', err_posterior); 
    end

    [valve.anterior pass err_anterior] = solve_valve_auto_continuation(valve.posterior, valve.tol_global, 'anterior'); 
    
    if pass 
        fprintf('Global solve passed anterior, err = %f\n\n', err_anterior); 
    else 
        fprintf('Global solve failed anterior, err = %f\n\n', err_anterior); 
    end 

    fig = figure; 
    fig = surf_plot(valve.posterior, fig);
    hold on 
    fig = surf_plot(valve.anterior, fig);
    title(sprintf('valve at p = %f', abs(p_0))); 
    
end
