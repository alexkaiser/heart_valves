function [leaflet pass err] = newton_solve_valve(leaflet, tol, max_it, max_consecutive_fails, max_total_fails) 
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

% some versions use an energy in addition to difference equations 
use_energy = false; 


% Get function handles for energy, diff eqns and Jacobian 
% if use_energy
%     energy   = leaflet.energy; 
% end 
diff_eqns = leaflet.diff_eqns; %@(leaflet) diff_eqns_and_linearize(leaflet, leaflet.diff_eqns); 
jacobian  = leaflet.jacobian; 


% Checks for a monotonic decrease if true 
% and decreases step length adaptively if not 
back_tracking = true; 
max_back_tracking_it = 20; 

% call 1D optimiazation routine if true 
optimization = false; 



consecutive_fails = 0; 
total_fails = 0; 

if ~exist('max_consecutive_fails', 'var')
    max_consecutive_fails = inf; 
end 

if ~exist('max_total_fails', 'var')
    max_total_fails = inf; 
end




plots = false; 
if plots 
    plot_freq = 1; 
    fig = figure; 
   
    surf_plot(leaflet, fig); 
    view(74,6); 
    hold off;  
end 

% newton step loop 
while err > tol
    
    tic 

    J = jacobian(leaflet); 

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
           
    
    F_linearized = diff_eqns(leaflet); 
    X_linearized_prev = linearize_internal_points(leaflet); 
    
    
    err_prev = err; 
    
    if use_energy 
        E_prev = E; 
    end 
    
    % solve the system,
    soln = J \ (-F_linearized); 

    % add in to get the next iterate 
    X_linearized = X_linearized_prev + soln; 

    % copy data back to 2d 
    leaflet = internal_points_to_2d(X_linearized, leaflet); 
    
    err = total_global_err(leaflet);         
    
    
    if back_tracking && (~optimization)
        
        alpha = 1.0; 
        back_tracking_it = 0; 
        while true 
           
            % pass
            if (err < err_prev) 
                fprintf('Line search passed with alpha = %f, \t ||F|| = %e\n', alpha, err); 
                consecutive_fails = 0; 
                break; 
            end 
            
            alpha = alpha / 2.0; 
            
            % add in to get the next iterate 
            X_linearized = X_linearized_prev + alpha * soln; 

            % copy data back to 2d 
            leaflet = internal_points_to_2d(X_linearized, leaflet); 

            err = total_global_err(leaflet); 
        
            back_tracking_it = back_tracking_it + 1; 
            
            if back_tracking_it > max_back_tracking_it
                warning('Failed to find a descent guess in allowed number of iterations.'); 
                consecutive_fails = consecutive_fails + 1; 
                total_fails = total_fails + 1; 
                
                % After search, if error has grown by an order or magnitude, return control
                if err > (10 * err_prev)
                    error('Error growth by over one order of magnitude. Return control to parent.'); 
                end 
                
                % Move along with this guess 
                break; 
            end 
        end 
        
    elseif back_tracking && optimization
        
        % norm squared of difference eqns 
        objective_fn = @(alpha) sum(diff_eqns(internal_points_to_2d(X_linearized_prev + alpha * soln, leaflet)).^2); 
        
        min_alpha = 0.0; 
        max_alpha = 1.0; 
        
        options = optimset('TolX', 1e2*tol, 'Display', 'off'); 
        
        [alpha_opt, norm_F_squared, exitflag, output] = fminbnd(objective_fn, min_alpha, max_alpha, options); 
        
        if exitflag == 1
            fprintf('One D optimization passed with alpha = %f, \t ||F|| = %e\n', alpha_opt, sqrt(norm_F_squared)); 
        else
            error('One dimensional optimization failed. Return control to parent.'); 
        end 
        
        X_linearized = X_linearized_prev + alpha_opt * soln; 
        leaflet = internal_points_to_2d(X_linearized, leaflet); 

    end 
    
    
    
%     
%     if back_tracking && use_energy
%         
%         % first order term in Taylor series, no coefficients
%         F_J_inv_F = F_linearized' * soln;  
%         
%         if F_J_inv_F <= 0.0
%            warning('Hessian has not made positive definite quadratic form');  
%         end 
%         
%         E = energy(leaflet); 
%         
%         alpha = 1.0; 
%         back_tracking_it = 0; 
%         while (E >= E_prev - c_backtrack * alpha * F_J_inv_F) 
%            
%             alpha = alpha / 2.0; 
%             
%             % add in to get the next iterate 
%             X_linearized = X_linearized_prev + alpha * soln; 
% 
%             % copy data back to 2d 
%             leaflet = internal_points_to_2d(X_linearized, leaflet); 
% 
%             E = energy(leaflet); 
%         
%             back_tracking_it = back_tracking_it + 1; 
%             
%             if back_tracking_it > max_back_tracking_it
%                 warning('failed to find a descent guess in allowed number of iterations'); 
%                 break; 
%             end 
%         end 
%         
%         err = total_global_err(leaflet); 
%     end 
%     
%     
%     if line_search && use_energy
%         
%         energy_gradient_hessian_handle = @(alpha) energy_gradient_hessian(X_linearized_prev + alpha*soln, leaflet); 
%         
%         min_alpha = 0.0; 
%         max_alpha = 1.0; 
%         
%         options = optimset('TolX', 1e2*tol, 'Display', 'off'); 
%         
%         [alpha_opt, E, exitflag, output] = fminbnd(energy_gradient_hessian_handle, min_alpha, max_alpha, options); 
%         
%         if exitflag == 1
%             fprintf('One D optimization passed with alpha = %e, \t E = %f\n', alpha_opt, E); 
%         else
%             error('One dimensional optimization failed'); 
%         end 
%         
%         X_linearized = X_linearized_prev + alpha_opt * soln; 
%         leaflet = internal_points_to_2d(X_linearized, leaflet); 
%     end 
%     
%     
    
    
    it = it + 1; 
    if it > max_it
        warning('Global solve failed to converge in %d iterations\n', it);
        pass = false; 
        break; 
    end  
    
    if isnan(err)
        error('NaN failure in solve. Return control to parent.'); 
    end 
    
    if consecutive_fails > max_consecutive_fails
        error('Too many consecutive failures on line search. Return control to parent.'); 
    end 
    
    if total_fails > max_total_fails
        error('Too many total failures on line search. Return control to parent.'); 
    end 
    
    
    if use_energy
        fprintf('Global iteration = %d, \tnorm %e, \tE = %e, \telapsed = %f\n', it, err, E, toc)
    else
        fprintf('Global iteration = %d, \tnorm %e, \telapsed = %f\n', it, err, toc)
    end 
    
    if plots && mod(it, plot_freq) == 0 
        surf_plot(leaflet, fig); 
        view(74,6); 
        title(sprintf('it = %d', it));         
        hold off; 
        pause(.1); 
    end 
    
    
    if isfield(leaflet, 'iteration_movie') && leaflet.iteration_movie 
        
        if ~isfield(leaflet, 'frame')
            error('must have frame number to make movie')
        end 
        
        
        output_leaflet_to_xyz_format(leaflet, ~leaflet.springs_written)
        
        leaflet.springs_written = true; 
        leaflet.frame = leaflet.frame + 1; 
        
    end    
end 

