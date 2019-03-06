function [x y] = svg_curve_to_points(curve_spec, x_0, y_0, min_x, max_x, min_y, max_y)
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

x_RELATIVE = zeros(length(curve_spec) / 6, 1); 
y_RELATIVE = zeros(length(curve_spec) / 6, 1); 

if mod(length(curve_spec),6) ~= 0
    error('must give 6 points per coordinate')
end

% points 5 6 are the accual coordinates 
idx = 1;
for i=5:6:length(curve_spec)
    x_RELATIVE(idx) = curve_spec(i); 
    y_RELATIVE(idx) = curve_spec(i+1); 
    idx = idx + 1; 
end 

x = zeros(length(x_RELATIVE)+1, 1); 
y = zeros(length(y_RELATIVE)+1, 1); 


x(1) = x_0; 
y(1) = y_0; 

for i=1:length(x_RELATIVE)
    x(i+1) = x(i) + x_RELATIVE(i); 
    y(i+1) = y(i) + y_RELATIVE(i); 
end 

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











