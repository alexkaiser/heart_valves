function [params pass err_over_time it] = solve_commissure_leaflet(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time, left)
%
% Commissural leaflet build 
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
%     left                  Left or right commissural leaflet 
% 
% Output 
%     params                Current solution parameters 
%     pass                  Whether solution has been computed to desired tolerance 
% 


pass = true; 

err = total_global_err_commissure(params, filter_params, left); 

it = start_it; 
fig = figure; 


while err > tol_global
    
    tic 
   
    % newton step here 

    % build the jacobian 
    J = sparse(build_jacobian_commissure(params, filter_params, left)); 

    F = difference_equations_commissure(params, filter_params, left); 
    F_linearized = linearize_internal_points_commissure(F, params); 
    X_linearized_prev = linearize_internal_points_commissure(params.X, params); 

    % solve the system,
    soln = J \ (-F_linearized);

    % add in to get the next iterate 
    X_linearized = X_linearized_prev + soln; 

    % copy data back to 2d 
    params = internal_points_to_2d_commissure(X_linearized, params); 
    % newton_step_coeff
    err = total_global_err_commissure(params, filter_params, left); 
    
    it = it + 1; 
    if it > max_it_global
        warning('Global solve failed to converge in %d iterations\n', it);
        pass = false; 
        break; 
    end  
    
    err_over_time(it) = err; 
    fprintf('Global iteration = %d, \tnorm %e, \telapsed = %f\n', it, err, toc)
        
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
        
        % input('stopped for fun...'); 
    end         
end 


if pass 
    data_name = sprintf('final_step_data_iteration_%d', it); 
    save(data_name); 
end 


































