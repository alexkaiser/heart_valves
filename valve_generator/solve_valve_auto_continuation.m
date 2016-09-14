function [leaflet pass err] = solve_valve_auto_continuation(leaflet, tol, max_it, name)
% 
% Automatically runs a continutation sequence 
% First tries the current reference fraction 
% If this fails then it is adaptively reduced and re-run
% 


[leaflet_current pass_current err] = newton_solve_valve(leaflet, tol, max_it);  

% quick exit if things work 
if pass_current 
    fprintf('Initial solve passed\n'); 
    leaflet = leaflet_current; 
    pass    = pass_current; 
    return
end 

fprintf('Initial solve failed, applying adaptive continuation\n'); 

% copy the last correct parameters 
leaflet_okay = leaflet; 

ref_current   = ref_frac / 2; 
ref_increment = ref_frac / 4; 

ever_passed = false; 


while true 
    
    fprintf('Solving with reference frac = %f\n', ref_current); 
    leaflet_okay.ref_frac = ref_current; 
    
    [leaflet_current pass_current err] = newton_solve_valve(leaflet_okay, tol, max_it);  
    
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
        leaflet_okay = leaflet_current; 
        
    else
        
        if ever_passed
            % if the current setup has passed, increment less  
            ref_increment = ref_increment / 4;  
            ref_current   = ref_last_passed + ref_increment; 
        else 
            % if never passed, current and incement must decrease 
            ref_increment = ref_increment / 4;  
            ref_current   = ref_current / 2; 
        end
        
    end 
end 
    
  
if pass_current 
    leaflet = leaflet_current; 
    pass    = pass_current; 
    return
else
    error('Auto continuation failed.'); 
end 












