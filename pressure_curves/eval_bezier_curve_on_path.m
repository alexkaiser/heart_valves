function [x y] = eval_bezier_curve_on_path(N_per_unit_length, curve_spec, x_0, y_0, min_x, max_x, min_y, max_y)
% 
% Takes coordinates from an svg curve and returns x,y pairs 
%  
% The 5th and 6th (and up in strides of 6) have the coordinates 
% 
% Initial coordinates must be provided 
% Movement is relative 
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

if mod(length(curve_spec),6) ~= 0
    error('must give 6 points per coordinate')
end

% The bezier curve is constructed with four points 

p0 = [x_0, y_0]; 

vals = zeros(length(curve_spec)/6 * N_per_unit_length,2); 

idx = 1; 

dt = 1/N_per_unit_length; 


for i = 1:6:length(curve_spec)

    % the positions are detmined relative to the current p0
    p1 = p0 + curve_spec(i  :i+1); 
    p2 = p0 + curve_spec(i+2:i+3); 
    p3 = p0 + curve_spec(i+4:i+5); 

    curve = @(t)    (1-t).^3          * p0 ...
                + 3*(1-t).^2 .* t     * p1 ...
                + 3*(1-t)    .* t.^2  * p2 ...
                +               t.^3  * p3; 
    
    % evaluate on mesh 
    vals(idx : idx+N_per_unit_length-1,:) = curve( (dt:dt:1)' );         
            
    idx = idx + N_per_unit_length; 
    
    p0 = p3; 
end 

% crop unused from ceiling function or whatever sloppy allocation 
vals = vals(1:idx-1,:); 

% just append the zero-th point on at the end here 
x = [x_0; vals(:,1)]; 
y = [y_0; vals(:,2)]; 

% the y coordinate is origin up, reverse it
y = -y; 


% zero the minimum, scale the array, fix the real minimum
x = x - min(x); 
range = max(x); 
range_actual = max_x - min_x; 
x = (range_actual / range) * x;
x = x + min_x; 

y = y - min(y); 
range = max(y); 
range_actual = max_y - min_y; 
y = (range_actual / range) * y; 
y = y + min_y; 









