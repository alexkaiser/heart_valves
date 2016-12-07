function [leaflet pass err] = newton_solve_valve(leaflet, tol, max_it) 
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
err = total_global_err(leaflet); 
it = 0; 

% Checks for a monotonic decrease if true 
% and decreases step length adaptively if not 
back_tracking = true; 
max_back_tracking_it = 50; 

% newton step loop 
while err > tol
    
    tic 

    J = leaflet.jacobian(leaflet); 

    jacobian_cond_info = false; 
    if jacobian_cond_info
        condition_num = cond(J)
    end 

    jacobian_det_info = false;  
    if jacobian_det_info
        'full determinant'
        det(J)

        'on diagonal blocks'
        total_internal = 3*sum(leaflet.is_internal(:)); 
        for block_start = 1:3:(total_internal - 1)
            
            col_rank = rank(full(J(:,block_start + (0:2))), tol);  
            
            if col_rank < 3
                'found a rank deficient 3x3 block of columns'
                block_start
                col_rank
            end 
            
        end 
    end 
           
    [F F_chordae_left F_chordae_right] = leaflet.diff_eqns(leaflet); 
    F_linearized      = linearize_internal_points(leaflet, F, F_chordae_left, F_chordae_right); 
    X_linearized_prev = linearize_internal_points(leaflet, leaflet.X, leaflet.chordae.C_left, leaflet.chordae.C_right); 
    
    
    err_prev = err; 
    
    % solve the system,
    soln = J \ (-F_linearized); 

    % add in to get the next iterate 
    X_linearized = X_linearized_prev + soln; 

    % copy data back to 2d 
    leaflet = internal_points_to_2d(X_linearized, leaflet); 
    
    err = total_global_err(leaflet);         
    
    
    if back_tracking 
        
        alpha = 1.0; 
        back_tracking_it = 0; 
        while (err > err_prev) 
           
            alpha = alpha / 2.0; 
            
            % add in to get the next iterate 
            X_linearized = X_linearized_prev + alpha * soln; 

            % copy data back to 2d 
            leaflet = internal_points_to_2d(X_linearized, leaflet); 

            err = total_global_err(leaflet); 
        
            back_tracking_it = back_tracking_it + 1; 
            
            if back_tracking_it > max_back_tracking_it
                warning('failed to find a decent guess in allowed number of iterations'); 
                break; 
            end 
        end 
        
    end 
    
    
    
    
    it = it + 1; 
    if it > max_it
        warning('Global solve failed to converge in %d iterations\n', it);
        pass = false; 
        break; 
    end  
     
    fprintf('Global iteration = %d, \tnorm %e, \telapsed = %f\n', it, err, toc)
    
    
    if isfield(leaflet, 'iteration_movie') && leaflet.iteration_movie 
        
        if ~isfield(leaflet, 'frame')
            error('must have frame number to make movie')
        end 
        
        
        output_leaflet_to_xyz_format(leaflet, ~leaflet.springs_written)
        
        leaflet.springs_written = true; 
        leaflet.frame = leaflet.frame + 1; 
        
    end    
end 

