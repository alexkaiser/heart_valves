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

% Copyright (c) 2019, Alexander D. Kaiser
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

pass = true; 
err = total_global_err(leaflet); 
it = 0; 

fprintf('Global iteration = %d, \tnorm %e\n', it, err)

% Get function handles for diff eqns and Jacobian 
diff_eqns = leaflet.diff_eqns; 
jacobian  = leaflet.jacobian; 


% Checks for a monotonic decrease if true 
% and decreases step length adaptively if not 
back_tracking = true; 
max_back_tracking_it = 20; 

% call 1D optimiazation routine if true 
optimization = false; 

consecutive_fails = 0; 
total_fails = 0; 

% will run a single iteration 
% residual computed by converting previous solution to sybolic datatype
% generally ineffective 
iterative_refinement = false; 

% if true, will continue to run iterations of Newton's method 
% until the error fails to decrease 
extra_iterations = true; 

% maximum number of extra iterations 
max_extra_iterations = 10; 
n_extra_it = 0; 

% solves linear system with optional advanpix package 
advanpix_multiprecision = false; 

% summary of current iteration information is printed as a format that goes
% into a latex table if true
% human-friendly format if false 
latex_table = false; 


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
    title('Initial');  
    view(74,6); 
    hold off;  
end 


if isfield(leaflet, 'iteration_movie') && leaflet.iteration_movie 
    
    movie_name = leaflet.movie_name; 
    
    fig_movie = figure; 
    surf_plot(leaflet, fig_movie); 
    grid off 
    axis off 
    set(fig_movie, 'Position', [100, 100, 1000, 1000])
    set(fig_movie,'PaperPositionMode','auto')
    title_str = sprintf('%s_diagonal_%d', movie_name, it); 
    printfig(fig_movie, title_str); 
    view(90,0)
    title_str = sprintf('%s_front_%d', movie_name, it); 
    printfig(fig_movie, title_str); 
    view(0,90)
    title_str = sprintf('%s_top_%d', movie_name, it);
    printfig(fig_movie, title_str); 
    % close(fig_movie); 
end 

converged = false; 
error_out = false;

% newton step loop 
while true
    
    if err < tol 
        if ~converged
            fprintf(' -- Newton: Solve converged to tolerance. Running extra iterations until error stops decreasing\n'); 
        end 
        
        converged = true; 
        leaflet_on_loop_entry = leaflet; 
        pass_on_loop_entry = pass; 
        err_on_loop_entry = err; 
        
        n_extra_it = n_extra_it + 1; 
        
        if n_extra_it > max_extra_iterations
            fprintf(' -- Newton: Max extra iterations reached'); 
            error_out = true; 
            break; 
        end 
        
    end 
    
    if converged && (~extra_iterations)
        break
    end 
    
    tic 

    J = jacobian(leaflet); 

    jacobian_cond_info = true; 
    if jacobian_cond_info
        condition_num = condest(J); 
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
 

    
    if advanpix_multiprecision 
        digits = 34; 
        mp.Digits(digits);
        
        J_extended_prec = mp(J); 
        F_linearized_extended_prec = mp(F_linearized); 
        soln_extended_prec = J_extended_prec \ (-F_linearized_extended_prec); 
        
%         err_abs = norm(soln_extended_prec - soln); 
%         err_rel = norm(soln_extended_prec - soln) / norm(soln_extended_prec); 
%         
%         fprintf('err_rel = %.20e\t', err_rel)
%         fprintf('err_abs = %.20e\n', err_abs)
        
        use_soln_cast_to_double = true; 
        if use_soln_cast_to_double
            soln = double(soln_extended_prec); 
        else
            fprintf(' -- Newton: NO multiprecision used\n')
        end
        
    else 
        % default precision 
        % solve the system,
        soln = J \ (-F_linearized); 
    end 
    
    
    if iterative_refinement
        % see
        % https://blogs.mathworks.com/cleve/2015/02/16/iterative-refinement-for-solutions-to-linear-systems/ 
        
        % simple triple precision residual 
        % residual = residual3p(J,soln,-F_linearized); 
        % redidual with double to sym to double 
        residual = double(J*sym(soln,'f') - (-F_linearized)); 
        correction = J \ residual; 
        
        soln = soln - correction; 
        
        check_iterative_refinement(J)
    end 
    
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
                % convergence on basic Newton step, no line search 
                % fprintf('Line search passed with alpha = %f, \t ||F|| = %e\n', alpha, err); 
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
                fprintf(' -- Warning:  Failed to find a descent guess in allowed number of iterations.\n');                 
                fprintf(' -- Newton: Difference in no descent guess with previous guess = %.20e\n', norm(X_linearized - X_linearized_prev)); 
                
                consecutive_fails = consecutive_fails + 1; 
                total_fails = total_fails + 1; 
                
                % After search, if error has grown by an order or magnitude, return control
                if err > (10 * err_prev)
                    error_out = true; 
                    warning('Error growth by over one order of magnitude. Return control to parent.'); 
                    break; 
                end 
                
                % Move along with this guess 
                break; 
            end 
            
            % double break if error_out has been set to true 
            if error_out 
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
        
    % various ways to leave the main loop 
    it = it + 1; 
    if it > max_it
        error_out = true; 
        fprintf(' -- Warning: Global solve failed to converge in %d iterations. Return control to parent.\n', it);
        break; 
    end  
    
    if isnan(err)
        error_out = true; 
        fprintf(' -- Warning: NaN failure in solve. Return control to parent.\n'); 
        break; 
    end 
    
    if consecutive_fails > max_consecutive_fails
        error_out = true; 
        fprintf(' -- Warning: Too many consecutive failures on line search. Return control to parent.\n'); 
        break; 
    end 
    
    if total_fails > max_total_fails
        error_out = true; 
        fprintf(' -- Warning: Too many total failures on line search. Return control to parent.\n'); 
        break; 
    end 
    

    % Summary of current iteration information 
    if latex_table
        fprintf('%d & %.2e & %f & %.2f \\\\ \n \\hline \n', it, err, alpha, toc); 
    else
        fprintf('Global iteration = %d, \tnorm %e, \telapsed = %f', it, err, toc); 
        if exist('alpha', 'var')
            fprintf(', \talpha = %f', alpha); 
        end
        if jacobian_cond_info 
            fprintf(', \tcond = %.1e', condition_num); 
        end 
        fprintf('\n'); 
    end 
        
    
    if plots && mod(it, plot_freq) == 0 
        fig = figure; 
        surf_plot(leaflet, fig); 
        view(74,6); 
        title(sprintf('it = %d', it));         
        hold off; 
        pause(.1); 
    end 
    
    
    if isfield(leaflet, 'iteration_movie') && leaflet.iteration_movie 
        
%         output_leaflet_to_xyz_format(leaflet, ~leaflet.springs_written)
%         
%         leaflet.springs_written = true; 

        fig_movie = figure; 
        surf_plot(leaflet, fig_movie); 
        grid off 
        axis off 
        set(fig_movie, 'Position', [100, 100, 1000, 1000])
        set(fig_movie,'PaperPositionMode','auto')
        title_str = sprintf('%s_diagonal_%d', movie_name, it); 
        printfig(fig_movie, title_str); 
        view(90,0)
        title_str = sprintf('%s_front_%d', movie_name, it); 
        printfig(fig_movie, title_str); 
        view(0,90)
        title_str = sprintf('%s_top_%d', movie_name, it);
        printfig(fig_movie, title_str); 
%         close(fig_movie); 
    end    
end 

if converged && ~(error_out)
    % all is well, no cleanup needed
    return; 

elseif error_out && converged 
    % left after extra iterations
    % solver converged then took a bad step
    % set previous good step to return values 
    
    fprintf(' -- Newton: In final step, err_rejected = %.20e, err_on_loop_entry = %.20e\n', err, err_on_loop_entry); 
    fprintf(' -- Newton: In final step, err_differences = %.20e\n', err - err_on_loop_entry); 
    
    fprintf(' -- Newton: Before taking old entry, ')
    if all(all(all( (leaflet_on_loop_entry.X == leaflet.X) | isnan(leaflet.X) )))
        fprintf('all equal in leaflets or nan masked\n')
    else 
        fprintf('NOT all equal in leaflets\n')
    end 
    
    leaflet = leaflet_on_loop_entry; 
    pass = pass_on_loop_entry; 
    err = err_on_loop_entry; 
    
    fprintf(' -- Newton: After taking old entry, ')
    if all(all(all( (leaflet_on_loop_entry.X == leaflet.X) | isnan(leaflet.X) )))
        fprintf('all equal or nan masked\n') 
    else 
        fpritnf('NOT all equal in leaflets\n')
    end 
    
elseif error_out && (~converged)
    % this is an actual failure for nonconvergence 
    error('Solver failed to converge');     
else     
    error('Should be impossible, check control flow'); 
end 



fprintf('Final iteration       \tnorm %e, \telapsed = %f', err, toc); 
if jacobian_cond_info 
    fprintf(', \tcond = %.1e', condition_num); 
end 
fprintf('\n'); 




