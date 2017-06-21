function [leaflet pass err any_passed] = solve_valve_pressure_auto_continuation(leaflet, tol, max_it, max_continuations, p_easy, p_goal, max_consecutive_fails, max_total_fails)
% 
% Automatically runs a continutation sequence 
% 
% If any solve has passed then the leaflet that passed closest to goal is returned 
% 
% 

initial_p_plot = true; 
plots = true; 


leaflet.p_0 = p_goal; 

try
    [leaflet_current pass err] = newton_solve_valve(leaflet, tol, max_it, max_consecutive_fails, max_total_fails);  

    % quick exit if fail 
    if pass 
        fprintf('Initial solve at goal passed.\n'); 
        leaflet = leaflet_current; 
        any_passed = true; 
        return; 
    else
        fprintf('Initial solve at goal failed. Running continuations.\n'); 
    end 
catch 
    fprintf('Initial solve at goal failed. Running continuations.\n'); 
end 


leaflet.p_0 = p_easy; 

[leaflet_current pass err] = newton_solve_valve(leaflet, tol, max_it, max_consecutive_fails, max_total_fails);  

% quick exit if fail 
if pass 
    fprintf('Solve at p initial passed. Continutations proceeding to goal.\n'); 
else
    error('Initial solve failed, no continuation possible. Adjust range.'); 
end 

fprintf('Applying adaptive continuation on pressure.\n\n'); 
any_passed = true; 

if initial_p_plot 
    fig = figure; 
    surf_plot(leaflet_current, fig); 
    title('Initial converged valve'); 
    pause(.1); 
end 

% copy the last correct parameters 
leaflet_okay = leaflet_current; 

p_last_passed = p_easy; 

p_current   = p_goal; 
p_increment = (p_goal - p_easy); 

if (p_goal - p_easy) < 0.0 
    increasing = false; 
else 
    increasing = true; 
end 

continuation_num = 1; 

while true 
    
    fprintf('Solving with p = %f\n', p_current); 
    
    leaflet_okay.p_0 = p_current; 
    pass = false; 
    
    try    
        [leaflet_current pass err] = newton_solve_valve(leaflet_okay, tol, max_it, max_consecutive_fails, max_total_fails);  
    
        if pass
        
            fprintf('Solve passed.\n\n', p_current); 
            
            if plots
                fig = figure; 
                surf_plot(leaflet_current, fig); 
                title(sprintf('Valve at p = %f', p_current)); 
                pause(.1); 
            end

            if p_current == p_goal
                leaflet_okay = leaflet_current; 
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
            continuation_num = continuation_num + 1; 
        end 
        
        
    catch 
        fprintf('\n\n'); 
        err = lasterror; 
        disp(err.message); 
        fprintf('Solve failed due to line search errors.\n')
        p_increment = p_increment / 4;  
        p_current   = p_last_passed + p_increment; 
        continuation_num = continuation_num + 1; 

    end 
    
    if continuation_num > max_continuations
        break
    end 
    
end 
    
leaflet = leaflet_okay; 
leaflet.p_0 = p_current; 

if ~pass 
    warning(sprintf('Leaving with pressure not equal to desired goal.\nFinal pressure = %f\n', leaflet_okay.p_0)); 
end 



