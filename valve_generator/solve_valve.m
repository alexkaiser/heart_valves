function [valve pass_all] = solve_valve(valve, p_range, repulsive_coeff_range)
% 
% Refines valve data structure to equilibrium 
% Applies auto-continuation to ref_frac and updates both leaflets 
% 


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
    
end 