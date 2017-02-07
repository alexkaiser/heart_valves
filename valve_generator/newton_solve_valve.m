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

% some versions use an energy in addition to difference equations 
use_energy = false; 


% Get function handles for energy, diff eqns and Jacobian 
if use_energy
    energy   = leaflet.energy; 
end 
diff_eqns = @(leaflet) diff_eqns_and_linearize(leaflet, leaflet.diff_eqns); 
jacobian  = leaflet.jacobian; 

if isfield(leaflet, 'repulsive_potential') && leaflet.repulsive_potential
    
    if use_energy 
        energy = @(leaflet) (energy(leaflet) + energy_repulsive(leaflet)); 
    end 
    diff_eqns  = @(leaflet) (diff_eqns(leaflet) + diff_eqns_and_linearize(leaflet, @difference_equations_repulsive)); 
    jacobian   = @(leaflet) (jacobian(leaflet)  + build_jacobian_replusive(leaflet)); 
end 


% Checks for a monotonic decrease if true 
% and decreases step length adaptively if not 
back_tracking = true; 
max_back_tracking_it = 12; 
if back_tracking 
    if use_energy
        use_energy = true; 
        E = energy(leaflet); 
        c_backtrack = .9; 
    end 
end 

line_search = false; 
if line_search
    if isfield(leaflet, 'energy')
        use_energy = true; 
        E = energy(leaflet); 
    else 
        error('Cannot line search without energy'); 
    end         
end 


if back_tracking && line_search
    error('only one line search strategy allowed'); 
end 


plots = true; 
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
    X_linearized_prev = linearize_internal_points(leaflet, leaflet.X, leaflet.chordae.C_left, leaflet.chordae.C_right); 
    
    
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
    
    
    if back_tracking && (~use_energy)
        
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
                warning('failed to find a descent guess in allowed number of iterations'); 
                break; 
            end 
        end 
        
    end 
    
    
    if back_tracking && use_energy
        
        % first order term in Taylor series, no coefficients
        F_J_inv_F = F_linearized' * soln;  
        
        if F_J_inv_F <= 0.0
           warning('Hessian has not made positive definite quadratic form');  
        end 
        
        E = energy(leaflet); 
        
        alpha = 1.0; 
        back_tracking_it = 0; 
        while (E >= E_prev - c_backtrack * alpha * F_J_inv_F) 
           
            alpha = alpha / 2.0; 
            
            % add in to get the next iterate 
            X_linearized = X_linearized_prev + alpha * soln; 

            % copy data back to 2d 
            leaflet = internal_points_to_2d(X_linearized, leaflet); 

            E = energy(leaflet); 
        
            back_tracking_it = back_tracking_it + 1; 
            
            if back_tracking_it > max_back_tracking_it
                warning('failed to find a descent guess in allowed number of iterations'); 
                break; 
            end 
        end 
        
        err = total_global_err(leaflet); 
    end 
    
    
    if line_search && use_energy
        
        energy_gradient_hessian_handle = @(alpha) energy_gradient_hessian(X_linearized_prev + alpha*soln, leaflet); 
        
        min_alpha = 0.0; 
        max_alpha = 1.0; 
        
        options = optimset('TolX', 1e2*tol, 'Display', 'off'); 
        
        [alpha_opt, E, exitflag, output] = fminbnd(energy_gradient_hessian_handle, min_alpha, max_alpha, options); 
        
        if exitflag == 1
            fprintf('One D optimization passed with alpha = %e, \t E = %f\n', alpha_opt, E); 
        else
            error('One dimensional optimization failed'); 
        end 
        
        X_linearized = X_linearized_prev + alpha_opt * soln; 
        leaflet = internal_points_to_2d(X_linearized, leaflet); 
    end 
    
    
    
    
    it = it + 1; 
    if it > max_it
        warning('Global solve failed to converge in %d iterations\n', it);
        pass = false; 
        break; 
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
        pause(0.5)
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

