function [params pass err_over_time] = solve_valve(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time)
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


err = total_global_err(); 

it = start_it; 
fig = figure; 

while err > tol_global
    
    tic 
   
    % newton step here 
    
    
    it = it + 1;
    
    if it > max_it_global
        warning('Global solve failed to converge in %d iterations\n', it);
        pass = false; 
        break; 
    end  
    
    err = total_global_err(  ); 
    
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
        
        if ~isempty(chordae)
            tree_plot(params, fig);
        end 
        
        %input('stopped for fun...'); 
            
    end 
        
end 


if pass 
    data_name = sprintf('final_step_data_iteration_%d', it); 
    save(data_name); 
end 


































