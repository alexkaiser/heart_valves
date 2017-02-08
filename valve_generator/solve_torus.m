function torus = solve_torus(torus, p_range)
% 
% Refines valve data structure to equilibrium 
% Applies auto-continuation to ref_frac and updates both leaflets 
% 



for p_0 = p_range

    % when the pressure changes, just update the pressure and re-run the setup 
    torus.p_0 = p_0; 

    [torus pass err] = solve_valve_auto_continuation(torus, torus.tol_global, torus.max_it, 'torus'); 

    if pass 
        fprintf('Global solve passed anterior, err = %f\n\n', err); 
    else 
        fprintf('Global solve failed anterior, err = %f\n\n', err); 
    end 

    fig = figure; 
    fig = torus_plot(torus, fig);
    title(sprintf('torus at p = %f', abs(p_0))); 
end 