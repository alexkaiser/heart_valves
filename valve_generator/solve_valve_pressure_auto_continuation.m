function [leaflet pass err] = solve_valve_pressure_auto_continuation(leaflet, tol, max_it_init, max_it_continuation, p_initial, p_goal, max_consecutive_fails, max_total_fails)
% 
% Automatically runs a continutation sequence 
% First tries the current reference fraction 
% If this fails then it is adaptively reduced and re-run
% 


leaflet.p_0 = p_initial; 

plots = true; 

[leaflet_current pass err] = newton_solve_valve(leaflet, tol, max_it_init, max_consecutive_fails, max_total_fails);  

% quick exit if fail 
if pass 
    fprintf('Initial solve passed.\n'); 
else
    error('Initial solve failed, no continuation possible. Adjust range.'); 
end 

fprintf('Applying adaptive continuation on pressure.\n\n'); 

if plots
    fig = figure; 
    surf_plot(leaflet_current, fig); 
    title('Initial converged valve'); 
    pause(.1); 
end 

% copy the last correct parameters 
leaflet_okay = leaflet_current; 

p_last_passed = p_initial; 

p_current   = p_goal; 
p_increment = (p_goal - p_initial); 

if (p_goal - p_initial) < 0.0 
    increasing = false; 
else 
    increasing = true; 
end 

while true
    
    fprintf('Solving with p = %f\n', p_current); 
    
    leaflet_okay.p_0 = p_current; 
    
    %try
        [leaflet_current pass err] = newton_solve_valve(leaflet_okay, tol, max_it_continuation, max_consecutive_fails, max_total_fails);  
    
        if pass
        
            fprintf('Solve passed.\n\n', p_current); 
            
            if plots
                fig = figure; 
                surf_plot(leaflet_current, fig); 
                title(sprintf('Valve at p = %f', p_current)); 
                pause(.1); 
            end

            if p_current == p_goal
                break
            end 

            p_last_passed = p_current; 

            % increment, but do not pass goal 

            if increasing 
                p_current = min(p_current + p_increment, p_goal);
            else 
                p_current = max(p_current + p_increment, p_goal);
            end 

            leaflet_okay = leaflet_current; 
        
        else 
            
            fprintf('Solve failed due to max iterations.\n')
            p_increment = p_increment / 4;  
            p_current   = p_last_passed + p_increment; 
        end 
        
        
%     catch 
%         fprintf('\n\n'); 
%         err = lasterror; 
%         disp(err.message); 
%         fprintf('Solve failed due to line search errors.\n')
%         p_increment = p_increment / 4;  
%         p_current   = p_last_passed + p_increment; 
% 
%     end 
end 
    

if pass 
    leaflet = leaflet_current; 
    return
else
    error('Auto continuation failed.'); 
end 



