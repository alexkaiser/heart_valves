function points = get_papillary_coords(valve, papillary_idx, n_points, min_angle, max_angle)
% 
% Returns a distributed list of papillary points 
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


papillary = valve.skeleton.papillary(papillary_idx); 
c         = papillary.center; 
r         = papillary.radius; 
n         = papillary.normal; 

% this points from left to right papillary
% this is approximately in the 'y' direction, but determined by tip location 
approx_vertical = valve.skeleton.l_to_r_papillary; 

% project onto orthogonal complement of span(n)
% and normalize 
% local 'y' coord in the circle 
approx_vertical = (eye(3) - n*n') * approx_vertical; 
approx_vertical = approx_vertical / norm(approx_vertical); 

% local 'x' coord in the circle 
approx_horiz    = cross(approx_vertical, n); 

if abs(norm(approx_horiz) - 1) > eps 
    error('Should have gotten a unit vector by construction');     
end 

angles = linspace(min_angle, max_angle, n_points); 
points = zeros(3,n_points); 

for i=1:n_points 
    points(:,i) = c + ...
                  r * cos(angles(i)) * approx_horiz + ...
                  r * sin(angles(i)) * approx_vertical; 
end 





