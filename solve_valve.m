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


pass = true; 

err = total_global_err(params, filter_params); 

it = start_it; 


while err > tol_global
    
    tic 
   

    % newton step here 

    % build the jacobian 
    J = sparse(build_jacobian(params, filter_params)); 

    jacobian_cond_info = false; 
    if jacobian_cond_info
        condition_num = cond(J)
    end 

    jacobian_det_info = false;  
    if jacobian_det_info
        'full determinant'
        det(J)

        'on diagonal blocks'
        total_internal = 3*params.N*(params.N+1)/2; 
        for block_start = 1:3:(total_internal - 1)

%                 if abs(det(J(block_start + (0:2), block_start + (0:2)))) < 1e-2
%                     block_start
%                     det(J(block_start + (0:2), block_start + (0:2)))
%                 end 

            tol = 1e2 * eps; 
            col_rank = rank(J(:,block_start + (0:2)), tol);  

            if col_rank < 3
                'found a rank deficient 3x3 block of columns'
                block_start
                col_rank
            end 


        end 

    end 

    F = difference_equations(params, filter_params); 
    F_linearized = linearize_internal_points(F, params); 
    X_linearized_prev = linearize_internal_points(params.X, params); 

    % solve the system,
    soln = J \ (-F_linearized);

    % add in to get the next iterate 
    X_linearized = X_linearized_prev + soln; 

    % copy data back to 2d 
    params = internal_points_to_2d(X_linearized, params); 

    err = total_global_err(params, filter_params);         
       
    
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
        if exist('fig', 'var') 
            if ishandle(fig)
                close(fig)
            end
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


































