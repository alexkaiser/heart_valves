function [params pass err_over_time it] = solve_valve_auto_continuation(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time, ref_frac, name)

% Automatically runs a continutation sequence 
% First tries the current reference fraction 
% If this fails then it is adaptively reduced and re-run

[params_current pass_current err_over_time it] = solve_valve(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time); 

% quick exit if things work 
if pass_current 
    fprintf('Initial solve passed\n'); 
    params = params_current; 
    pass   = pass_current; 
    return
end 

fprintf('Initial solve failed, applying adaptive continuation\n'); 

% copy the last correct parameters 
params_okay = params; 

ref_current   = ref_frac / 2; 
ref_increment = ref_frac / 4; 

ever_passed = false; 


while true 
    
    fprintf('Solving with reference frac = %f\n', ref_current); 
    
    params_okay.ref_frac = ref_current; 
    
    [params_current pass_current err_over_time it] = solve_valve(params_okay, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time); 
    
    if pass_current 
        
        fprintf('Solve passed\n\n'); 
        
        if exist('name', 'var')
            data_name = sprintf('%s_at_ref_frac_%f', name, ref_current); 
            save(data_name); 
        end
        
        if ref_current == ref_frac
            break
        end 
        
        ever_passed = true;
        ref_last_passed = ref_current; 
        
        % increment, but do not pass goal 
        ref_current = min(ref_current + ref_increment, ref_frac);
        params_okay = params_current; 
        
 
            
    else
        
        if ever_passed
            % if the current setup has passed, just incremet less  
            ref_increment = ref_increment / 4;  
            ref_current   = ref_last_passed + ref_increment; 
        else 
            ref_increment = ref_increment / 4;  
            ref_current   = ref_current / 2; 
        end
        
    end 
end 
    
  
if pass_current 
    params = params_current; 
    pass   = pass_current; 
    return
else
    error('Auto continuation failed.'); 
end 












