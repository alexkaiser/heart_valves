function [valve valve_linear pass_all] = solve_valve(valve, p_range, linear_open_config, p_range_linear, strain, repulsive_coeff_range)
% 
% Refines valve data structure to equilibrium 
% Applies auto-continuation to pressure and updates both leaflets 
% 

pass_all = true; 
fig = figure; 

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
p_goal    = valve.anterior.p_0; 

[valve.anterior pass_anterior err_anterior] = solve_valve_pressure_auto_continuation(valve.anterior, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

if pass_anterior 
    fprintf('Global solve passed anterior, err = %e\n\n', err_anterior); 
else 
    fprintf('Global solve failed anterior, err = %e\n\n', err_anterior); 
end 


surf_plot(valve.anterior, fig); 


p_initial = 0; 
p_goal    = valve.posterior.p_0; 

[valve.posterior pass_posterior err_posterior] = solve_valve_pressure_auto_continuation(valve.posterior, valve.tol_global, valve.max_it, valve.max_it_continuation, p_initial, p_goal, valve.max_consecutive_fails, valve.max_total_fails); 

if pass_anterior 
    fprintf('Global solve passed anterior, err = %e\n\n', err_posterior); 
else 
    fprintf('Global solve failed anterior, err = %e\n\n', err_posterior); 
end 

surf_plot(valve.posterior, fig); 

pass_all = pass_all && pass_anterior && pass_posterior; 






if pass_all

    fprintf('Closed configurations passed, generating open configuration with linear constitutive law.\n'); 
    valve_linear = valve; 
    valve_linear.anterior  = set_rest_lengths_and_constants_linear(valve.anterior,  strain); 
    valve_linear.posterior = set_rest_lengths_and_constants_linear(valve.posterior, strain);       
    valve_linear.diff_eqns = @difference_equations_linear; 
    valve_linear.jacobian  = @build_jacobian_linear; 


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
        
        valve_linear.comm_left  = set_rest_lengths_and_constants_linear(valve.comm_left,  strain); 
        valve_linear.comm_right = set_rest_lengths_and_constants_linear(valve.comm_right, strain);

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
    
else 
    error('Solves failed to converge, did not produce models.'); 
end 



return; 



debug = false; 

if isfield(valve, 'optimization') && valve.optimization 

    % less ambitious for now 
    valve.tol_global = 1e-5; 
    
    valve.anterior = solve_valve_optimization(valve.anterior, valve.tol_global); 
    
    fig = figure; 
    fig = surf_plot(valve.anterior, fig);
    
    valve.posterior = solve_valve_optimization(valve.posterior, valve.tol_global); 
    
    fig = figure; 
    fig = surf_plot(valve.posterior, fig);    
    
    'cool'; 
    
    pass_all = false; 
    
else 

    for p_0 = p_range

        % when the pressure changes, just update the pressure and re-run the setup 
        valve.anterior.p_0 = p_0; 
        valve.posterior.p_0 = p_0; 

        if valve.posterior.reflect_pressure
            valve.posterior.p_0 = -p_0; 
        end 
        
        [valve.anterior pass_anterior err_anterior] = solve_valve_auto_continuation(valve.anterior, valve.tol_global, valve.max_it, 'anterior'); 

        if pass_anterior 
            fprintf('Global solve passed anterior, err = %f\n\n', err_anterior); 
        else 
            fprintf('Global solve failed anterior, err = %f\n\n', err_anterior); 
        end 
        
        valve.anterior.tension_debug = true; 
        difference_equations_bead_slip(valve.anterior); 
        valve.anterior.tension_debug = false; 

        [valve.posterior pass_posterior err_posterior] = solve_valve_auto_continuation(valve.posterior, valve.tol_global, valve.max_it, 'posterior'); 
    
        if pass_posterior 
            fprintf('Global solve passed posterior, err = %f\n\n', err_posterior); 
        else 
            fprintf('Global solve failed posterior, err = %f\n\n', err_posterior); 
        end
        
        if isfield(valve, 'comm_left') && isfield(valve, 'comm_right')
            
            valve.comm_left.p_0 = p_0; 
            
            [valve.comm_left pass_comm_left err_comm_left] = solve_valve_auto_continuation(valve.comm_left, valve.tol_global, valve.max_it, 'comm_left');

            if pass_comm_left 
                fprintf('Global solve passed comm_left, err = %f\n\n', err_comm_left); 
            else 
                fprintf('Global solve failed comm_left, err = %f\n\n', err_comm_left); 
            end
        
            valve.comm_right.p_0 = p_0; 
            
            [valve.comm_right pass_comm_right err_comm_right] = solve_valve_auto_continuation(valve.comm_right, valve.tol_global, valve.max_it, 'comm_right');

            if pass_comm_right 
                fprintf('Global solve passed comm_right, err = %f\n\n', err_comm_right); 
            else 
                fprintf('Global solve failed comm_right, err = %f\n\n', err_comm_right); 
            end
            
        end 
        

        if debug
            valve.posterior.tension_debug = true; 
            difference_equations_bead_slip(valve.posterior); 
            valve.posterior.tension_debug = false; 
        end 
    
        fig = figure; 
        fig = valve_plot(valve, fig);
        title(sprintf('valve at p = %f', abs(p_0))); 

    end

    pass_all = pass_anterior && pass_posterior; 
    
    if pass_all 
        if exist('repulsive_coeff_range', 'var') 
            
            c_repulsive_circumferential_orig = valve.c_repulsive_circumferential; 
            c_repulsive_radial_orig          = valve.c_repulsive_radial; 
            c_repulsive_chordae_orig         = valve.c_repulsive_chordae; 
            
            
            for repulsive_coeff_multiplier = repulsive_coeff_range
            
                % update the three repulsive coefficients on both leaflets 
                valve.anterior.c_repulsive_circumferential = repulsive_coeff_multiplier * c_repulsive_circumferential_orig; 
                valve.anterior.c_repulsive_radial          = repulsive_coeff_multiplier * c_repulsive_radial_orig; 
                valve.anterior.c_repulsive_chordae         = repulsive_coeff_multiplier * c_repulsive_chordae_orig; 
                
                valve.posterior.c_repulsive_circumferential = repulsive_coeff_multiplier * c_repulsive_circumferential_orig; 
                valve.posterior.c_repulsive_radial          = repulsive_coeff_multiplier * c_repulsive_radial_orig; 
                valve.posterior.c_repulsive_chordae         = repulsive_coeff_multiplier * c_repulsive_chordae_orig; 
                
                                
                [valve.anterior pass_anterior err_anterior] = solve_valve_auto_continuation(valve.anterior, valve.tol_global, valve.max_it, 'anterior'); 

                if pass_anterior 
                    fprintf('Global solve passed anterior, err = %f\n\n', err_anterior); 
                else 
                    fprintf('Global solve failed anterior, err = %f\n\n', err_anterior); 
                end 

                [valve.posterior pass_posterior err_posterior] = solve_valve_auto_continuation(valve.posterior, valve.tol_global, valve.max_it, 'posterior'); 

                if pass_posterior 
                    fprintf('Global solve passed posterior, err = %f\n\n', err_posterior); 
                else 
                    fprintf('Global solve failed posterior, err = %f\n\n', err_posterior); 
                end

                fig = figure; 
                fig = surf_plot(valve.posterior, fig);
                hold on 
                fig = surf_plot(valve.anterior, fig);
                title(sprintf('valve at repulsive coeff scaling = %f', repulsive_coeff_multiplier)); 
        
            end 
        
        end 
    end 
    
    pass_all = pass_anterior && pass_posterior; 
    
    if pass_all
        
        fprintf('Closed configurations passed, generating open configuration with linear constitutive law.\n'); 
        
        valve_linear = valve; 
        valve_linear.anterior  = set_rest_lengths_and_constants_linear(valve.anterior,  strain, valve.left_papillary_diastolic, valve.right_papillary_diastolic); 
        valve_linear.posterior = set_rest_lengths_and_constants_linear(valve.posterior, strain, valve.left_papillary_diastolic, valve.right_papillary_diastolic);       
        valve_linear.diff_eqns = @difference_equations_linear; 
        valve_linear.jacobian  = @build_jacobian_linear; 
        
        if isfield(valve_linear, 'comm_left') && isfield(valve_linear, 'comm_right')
            valve_linear.comm_left  = set_rest_lengths_and_constants_linear(valve.comm_left,  strain, valve.left_papillary_diastolic,  valve.left_papillary_diastolic); 
            valve_linear.comm_right = set_rest_lengths_and_constants_linear(valve.comm_right, strain, valve.right_papillary_diastolic, valve.right_papillary_diastolic);
        end 
        
        if linear_open_config 

            for p_0 = p_range_linear

                % when the pressure changes, just update the pressure and re-run the setup 
                valve_linear.anterior.p_0 = p_0; 
                valve_linear.posterior.p_0 = p_0; 

                if valve_linear.posterior.reflect_pressure
                    valve_linear.posterior.p_0 = -p_0; 
                end 

                [valve_linear.anterior pass_anterior err_anterior] = solve_valve_auto_continuation(valve_linear.anterior, valve_linear.tol_global, valve_linear.max_it, 'anterior'); 

                if pass_anterior 
                    fprintf('Global solve passed anterior, err = %f\n\n', err_anterior); 
                else 
                    fprintf('Global solve failed anterior, err = %f\n\n', err_anterior); 
                end 

                [valve_linear.posterior pass_posterior err_posterior] = solve_valve_auto_continuation(valve_linear.posterior, valve_linear.tol_global, valve_linear.max_it, 'posterior'); 

                if pass_posterior 
                    fprintf('Global solve passed posterior, err = %f\n\n', err_posterior); 
                else 
                    fprintf('Global solve failed posterior, err = %f\n\n', err_posterior); 
                end

                valve_linear.posterior.tension_debug = true; 
                difference_equations_bead_slip(valve_linear.posterior); 
                valve_linear.posterior.tension_debug = false; 


                if isfield(valve_linear, 'comm_left') && isfield(valve_linear, 'comm_right')

                    valve_linear.comm_left.p_0 = p_0; 

                    [valve_linear.comm_left pass_comm_left err_comm_left] = solve_valve_auto_continuation(valve_linear.comm_left, valve_linear.tol_global, valve_linear.max_it, 'comm_left');

                    if pass_comm_left 
                        fprintf('Global solve passed comm_left, err = %f\n\n', err_comm_left); 
                    else 
                        fprintf('Global solve failed comm_left, err = %f\n\n', err_comm_left); 
                    end

                    valve_linear.comm_right.p_0 = p_0; 

                    [valve_linear.comm_right pass_comm_right err_comm_right] = solve_valve_auto_continuation(valve_linear.comm_right, valve_linear.tol_global, valve_linear.max_it, 'comm_right');

                    if pass_comm_right 
                        fprintf('Global solve passed comm_right, err = %f\n\n', err_comm_right); 
                    else 
                        fprintf('Global solve failed comm_right, err = %f\n\n', err_comm_right); 
                    end

                end 
                
                
                fig = figure; 
                fig = valve_plot(valve_linear, fig);
                title(sprintf('valve, linear constitutive, at p = %f', abs(p_0))); 

            end
        end 
    end
    
    pass_all = pass_anterior && pass_posterior;
    
end 


