function [leaflet pass err any_passed] = solve_valve_pressure_auto_continuation(leaflet, tol, max_it, max_continuations, p_easy, p_goal, max_consecutive_fails, max_total_fails, goal_first)
% 
% Automatically runs a continutation sequence 
% 
% If any solve has passed then the leaflet that passed closest to goal is returned 
% 
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


initial_p_plot = true; 
plots = true; 

if ~exist('goal_first', 'var')
    goal_first = true; 
end 

if goal_first
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

            % update pressure before break 
            p_last_passed = p_current; 
            
            if p_current == p_goal
                leaflet_okay = leaflet_current; 
                break
            end 

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
        error_obj = lasterror; 
        disp(error_obj.message); 
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
leaflet.p_0 = p_last_passed; 

if pass 
    if p_last_passed ~= p_current 
        error("solve passed but current pressure not equal to last passed pressure, control flow bug"); 
    end 
end 

if ~pass 
    warning(sprintf('Leaving with pressure not equal to desired goal.\nFinal pressure = %f\n', leaflet_okay.p_0)); 
end 



