function [params pass err_over_time it] = solve_valve(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time)
%
% Full valve build. 
% Solves the nonlinear difference equations at each component. 
% Then checks global L2 norm. 
% Repeats until below global tolerance 
% 
% Input 
%     params                Current solution parameters 
%     filter_params         Coffee cone filter parameters 
%     tol_global            Error tolerance for global solve 
%     max_it_global         Maximum global iterations 
%     plot_and_save_freq    Saves all data every this many iterations 
%     start_it              Initial iteration, may not be zero if the run is restarted  
%     err_over_time         Error through iterations 
% 
% Output 
%     params                Current solution parameters 
%     pass                  Whether solution has been computed to desired tolerance 
% 

full_newton = false; 


pass = true; 

err = total_global_err(params, filter_params); 

it = start_it; 
fig = figure; 

while err > tol_global
    
    tic 
   
    if full_newton 
        % newton step here 

        % build the jacobian 
        J = build_jacobian(params, filter_params); 

        jacobian_det_info = true 
        if jacobian_det_info
            'full determinant'
            det(J)

            'on diagonal blocks'
            total_internal = 3*params.N*(params.N+1)/2; 
            for block_start = 1:3:(total_internal - 1)

                if abs(det(J(block_start + (0:2), block_start + (0:2)))) < 1e-2
                    block_start
                    det(J(block_start + (0:2), block_start + (0:2)))
                end 

            end 

        end 

        F = difference_equations(params, filter_params); 
        F_linearized = linearize_internal_points(F, params); 

        if it == 0 
            spy(J)
            title('jacobian nonzero pattern on iteration zero');

            % open a new figure 
            fig = figure;
        end 

        % solve the system,
        soln = J \ (-F_linearized); 

        % add in to get the next iterate 
        F_linearized = F_linearized + soln; 
        
        % copy data back to 2d 
        params = internal_points_to_2d(F_linearized, params); 
        
    else 
        tol_local = 1e-3 * tol_global; 
        max_it_local = 100; 
        params = update_leaflet_red_black(params, filter_params, tol_local, max_it_local); 
    end 
    
    
    
    
    it = it + 1; 
    if it > max_it_global
        warning('Global solve failed to converge in %d iterations\n', it);
        pass = false; 
        break; 
    end  
    
    
    err = total_global_err(params, filter_params); 

    err_over_time(it) = err; 
    fprintf('Global iteration = %d, \tnorm %g, \telapsed = %f\n', it, err, toc)
        
    if mod(it, plot_and_save_freq) == 0
        
        data_name = sprintf('data_iteration_%d', it); 
        save(data_name); 
        
        % plot for entertainment... 
        
        % close if open
        if ishandle(fig)
            close(fig)
        end 
        
        fig = surf_plot(params, filter_params); 
        title_str = sprintf('Surface at iteration %d', it); 
        title(title_str); 
        
        input('stopped for fun...'); 
    end         
end 


if pass 
    data_name = sprintf('final_step_data_iteration_%d', it); 
    save(data_name); 
end 


































