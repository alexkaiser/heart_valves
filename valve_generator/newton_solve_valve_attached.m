function [valve pass err] = newton_solve_valve_attached(valve, tol, max_it) 
%
% Newton solve on attached 
% 

pass = true; 
err = total_global_err_attached(valve); 
it = 0; 

% Checks for a monotonic decrease if true 
% and decreases step length adaptively if not 
back_tracking = true; 
max_back_tracking_it = 16; 


plots = true; 
if plots 
    plot_freq = 10; 
    fig = figure; 
    valve_plot(valve, fig); 
    view(74,6); 
    pause(1); 
    hold off;  
end 


% newton step loop 
while err > tol
    
    tic 

    J = build_jacobian_bead_slip_attached(valve); 

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
           
    [F_anterior F_posterior F_chordae_left F_chordae_right] = difference_equations_bead_slip_attached(valve); 
    F_linearized      = linearize_internal_points_bead_slip_attached(valve, F_anterior, F_posterior, F_chordae_left, F_chordae_right); 
    X_linearized_prev = linearize_internal_points_bead_slip_attached(valve, valve.anterior.X, valve.posterior.X, valve.anterior.chordae.C_left, valve.anterior.chordae.C_right); 
    
    err_prev = err; 
    
    % solve the system,
    soln = J \ (-F_linearized); 

    % add in to get the next iterate 
    X_linearized = X_linearized_prev + soln; 

    % copy data back to 2d 
    valve = internal_points_to_2d_attached(X_linearized, valve); 
    
    err = total_global_err_attached(valve);
    
    if back_tracking 
        
        alpha = 1.0; 
        back_tracking_it = 0; 
        while (err > err_prev) 
           
            alpha = alpha / 2.0; 
            
            % add in to get the next iterate 
            X_linearized = X_linearized_prev + alpha * soln; 

            % copy data back to 2d 
            valve = internal_points_to_2d_attached(X_linearized, valve); 
    
            err = total_global_err_attached(valve); 
        
            back_tracking_it = back_tracking_it + 1; 
            
            if back_tracking_it > max_back_tracking_it
                warning('failed to find a descent guess in allowed number of iterations'); 
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
    
    if plots && mod(it, plot_freq) == 0 
        valve_plot(valve, fig); 
        view(74,6); 
        title(sprintf('it = %d', it));         
        hold off; 
        pause(1); 
    end 
    
%     if isfield(leaflet, 'iteration_movie') && leaflet.iteration_movie 
%         
%         if ~isfield(leaflet, 'frame')
%             error('must have frame number to make movie')
%         end 
%         
%         
%         output_leaflet_to_xyz_format(leaflet, ~leaflet.springs_written)
%         
%         leaflet.springs_written = true; 
%         leaflet.frame = leaflet.frame + 1; 
%         
%     end    
end 

