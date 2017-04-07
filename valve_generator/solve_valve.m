function [valve valve_linear pass_all] = solve_valve(valve, strain)
% 
% Refines valve data structure to equilibrium 
% Applies auto-continuation to pressure and updates both leaflets 
% 

pass_all = true; 
fig = figure; 


if isfield(valve, 'leaflet')
    
    p_initial = 0; 
    p_goal    = valve.leaflet.p_0; 

    [valve.leaflet pass_leaflet err_leaflet] = solve_valve_pressure_auto_continuation(valve.leaflet, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

    if pass_leaflet 
        fprintf('Global solve passed, err = %e\n\n', err_leaflet); 
    else 
        fprintf('Global solve failed, err = %e\n\n', err_leaflet); 
    end 

    surf_plot(valve.leaflet, fig); 
    pause(0.01);
    
else 


    if isfield(valve, 'comm_left') && isfield(valve, 'comm_right')

        p_initial = 0; 
        p_goal    = valve.comm_right.p_0; 

        [valve.comm_right pass_comm_right err_comm_right] = solve_valve_pressure_auto_continuation(valve.comm_right, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

        if pass_comm_right 
            fprintf('Global solve passed comm_right, err = %e\n\n', err_comm_right); 
        else 
            fprintf('Global solve failed comm_right, err = %e\n\n', err_comm_right); 
        end

        surf_plot(valve.comm_right, fig); 
        hold on

        p_initial = 0; 
        p_goal    = valve.comm_left.p_0;

        [valve.comm_left pass_comm_left err_comm_left] = solve_valve_pressure_auto_continuation(valve.comm_left, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

        if pass_comm_left 
            fprintf('Global solve passed comm_left, err = %e\n\n', err_comm_left); 
        else 
            fprintf('Global solve failed comm_left, err = %e\n\n', err_comm_left); 
        end

        surf_plot(valve.comm_left, fig); 

        pass_all = pass_all && pass_comm_left && pass_comm_right; 

    end 


    p_initial = 0; 
    p_goal    = valve.posterior.p_0; 

    [valve.posterior pass_posterior err_posterior] = solve_valve_pressure_auto_continuation(valve.posterior, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

    if pass_posterior 
        fprintf('Global solve passed posterior, err = %e\n\n', err_posterior); 
    else 
        fprintf('Global solve failed posterior, err = %e\n\n', err_posterior); 
    end 

    surf_plot(valve.posterior, fig); 
    pause(0.01); 


    p_initial = 0; 
    p_goal    = valve.anterior.p_0; 

    [valve.anterior pass_anterior err_anterior] = solve_valve_pressure_auto_continuation(valve.anterior, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

    if pass_anterior 
        fprintf('Global solve passed anterior, err = %e\n\n', err_anterior); 
    else 
        fprintf('Global solve failed anterior, err = %e\n\n', err_anterior); 
    end 


    surf_plot(valve.anterior, fig); 




    pass_all = pass_all && pass_anterior && pass_posterior; 

    valve_linear = []; 

end 

return; 



if pass_all

    fprintf('Closed configurations passed, generating open configuration with linear constitutive law.\n'); 
    valve_linear = valve; 
    valve_linear.anterior  = set_rest_lengths_and_constants(valve.anterior,  strain); 
    valve_linear.posterior = set_rest_lengths_and_constants(valve.posterior, strain);       
    valve_linear.diff_eqns = @difference_equations_linear; 
    valve_linear.jacobian  = @build_jacobian_linear; 

    if isfield(valve_linear, 'comm_left') && isfield(valve_linear, 'comm_right')    
        valve_linear.comm_left  = set_rest_lengths_and_constants(valve.comm_left,  strain); 
        valve_linear.comm_right = set_rest_lengths_and_constants(valve.comm_right, strain);
    end    

    
    if ~valve.collagen_constitutive
    
        p_initial = valve_linear.anterior.p_0; 
        p_goal    = 0; 

        [valve_linear.anterior pass_anterior err_anterior] = solve_valve_pressure_auto_continuation(valve_linear.anterior, valve_linear.tol_global, valve_linear.max_it, valve_linear.max_it_continuation, p_initial, p_goal, valve_linear.max_consecutive_fails, valve_linear.max_total_fails); 

        if pass_anterior 
            fprintf('Global solve passed anterior, err = %e\n\n', err_anterior); 
        else 
            fprintf('Global solve failed anterior, err = %e\n\n', err_anterior); 
        end 


        p_initial = valve_linear.posterior.p_0; 
        p_goal    = 0; 

        [valve_linear.posterior pass_posterior err_posterior] = solve_valve_pressure_auto_continuation(valve_linear.posterior, valve_linear.tol_global, valve_linear.max_it, valve_linear.max_it_continuation, p_initial, p_goal, valve_linear.max_consecutive_fails, valve_linear.max_total_fails); 

        if pass_anterior 
            fprintf('Global solve passed anterior, err = %e\n\n', err_posterior); 
        else 
            fprintf('Global solve failed anterior, err = %e\n\n', err_posterior); 
        end 


        if isfield(valve_linear, 'comm_left') && isfield(valve_linear, 'comm_right')

            p_initial = valve_linear.comm_left.p_0;
            p_goal    = 0;

            [valve_linear.comm_left pass_comm_left err_comm_left] = solve_valve_pressure_auto_continuation(valve_linear.comm_left, valve_linear.tol_global, valve_linear.max_it, valve_linear.max_it_continuation, p_initial, p_goal, valve_linear.max_consecutive_fails, valve_linear.max_total_fails); 

            if pass_comm_left 
                fprintf('Global solve passed comm_left, err = %f\n\n', err_comm_left); 
            else 
                fprintf('Global solve failed comm_left, err = %f\n\n', err_comm_left); 
            end

            p_initial = valve_linear.comm_left.p_0;
            p_goal    = 0; 

            [valve_linear.comm_right pass_comm_right err_comm_right] = solve_valve_pressure_auto_continuation(valve_linear.comm_right, valve_linear.tol_global, valve_linear.max_it, valve_linear.max_it_continuation, p_initial, p_goal, valve_linear.max_consecutive_fails, valve_linear.max_total_fails); 

            if pass_comm_right 
                fprintf('Global solve passed comm_right, err = %f\n\n', err_comm_right); 
            else 
                fprintf('Global solve failed comm_right, err = %f\n\n', err_comm_right); 
            end

            pass_all = pass_all && pass_comm_left && pass_comm_right; 

        end 
    end 
    
else 
    error('Solves failed to converge, did not produce models.'); 
end 

