function [leaflet pass err] = solve_valve_auto_continuation(leaflet, tol, max_it, name)
% 
% Automatically runs a continutation sequence 
% First tries the current reference fraction 
% If this fails then it is adaptively reduced and re-run
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

[leaflet_current pass_current err] = newton_solve_valve(leaflet, tol, max_it);  

% quick exit if things work 
if pass_current 
    fprintf('Initial solve passed\n'); 
    leaflet = leaflet_current; 
    pass    = pass_current; 
    return
else
    if isfield(leaflet, 'leaflet_only') && leaflet.leaflet_only 
        error('Initial solve failed on leaflet only version, no continuation possible'); 
    end 
end 

fprintf('Initial solve failed, applying adaptive continuation\n'); 

% copy the last correct parameters 
leaflet_okay = leaflet; 

if isfield(leaflet, 'frame')
    leaflet_okay.frame = leaflet_current.frame; 
end 

ref_current   = leaflet.ref_frac / 2; 
ref_increment = leaflet.ref_frac / 4; 

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
        
        if ref_current == leaflet.ref_frac
            break
        end 
        
        ever_passed = true;
        ref_last_passed = ref_current; 
        
        % increment, but do not pass goal 
        ref_current = min(ref_current + ref_increment, leaflet.ref_frac);
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
        
        if isfield(leaflet, 'frame')
            leaflet_okay.frame = leaflet_current.frame; 
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












