function skeleton = valve_points_ct_diastole(output)
% 
% Takes a center, radius and papillary coords 
% Moves center to origin 
% Takes papillary tips in rigid manner
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

if ~exist('output', 'var')
    output = false; 
end 


radius = 0.1 * 16.065878777687722; 

normal = [-0.70870559417192203; 0.61353537460948204; 0.34829689187850299]; 


% initial coordinates are in mm
% multiply by .1 everywhere to get cm 

center = 0.1 * [43.211190372043184; 12.739459457869332; -18.011383112933899]; 
papillary_left  = 0.1 * [62.051506; 1.486666; -45.702660]; 
papillary_right = 0.1 * [84.484848; 3.326258; -21.130402]; 


% translate center to origin 
papillary_left  = papillary_left  - center; 
papillary_right = papillary_right - center; 

% get angle to rotate normal to in the y,z plane 
% and rotate accordingly 
theta = atan2(normal(2), normal(1)); 
R = rotation_matrix_z(-theta); 

normal          = R*normal; 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 

if abs(normal(2)) > eps
    error('did not land in the xz plane'); 
end 

% get angle to place on positive z axis 
phi = atan2(normal(3), normal(1)); 
R = rotation_matrix_y(pi/2 - phi); 

normal          = R*normal; 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 

if (abs(normal(1)) > eps) || (abs(normal(2)) > eps) 
    error('did not land in the z axis'); 
end 

if abs(normal(3) - 1) > eps
    error('z component of end normal is not in the right place'); 
end 

normal_to_papillary_plane = cross(papillary_left, papillary_right); 

% compute the midpoint
% rotate this to the negative x axis 
midpoint = 0.5 * (papillary_left(1:2) + papillary_right(1:2)); 
theta = atan2(midpoint(2), midpoint(1)); 

R = rotation_matrix_z(-theta + pi); 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 


if output 

    radius 
    papillary_left
    papillary_right


    th = 0:.1:2*pi; 
    figure; 
    plot3(radius*cos(th), radius*sin(th), zeros(size(th))); 
    hold on 
    plot3(papillary_left(1), papillary_left(2), papillary_left(3), '*'); 
    plot3(papillary_right(1), papillary_right(2), papillary_right(3), 'o'); 
    legend('ring', 'left', 'right'); 
    title('mitral skeleton geometry'); 
    xlabel('x')
    ylabel('y')
    zlabel('z')

end 

l_to_r_papillary = (papillary_right - papillary_left); 
l_to_r_papillary = l_to_r_papillary / norm(l_to_r_papillary);


papillary_radius = 0; %0.25; 

warning('Diastolic papillary radius set to zero here.'); 

left_papillary_center  = papillary_left  + papillary_radius * l_to_r_papillary; 
right_papillary_center = papillary_right - papillary_radius * l_to_r_papillary; 



% always horizontal here 
normal = [0; 0; 1]; 

skeleton.r                   = radius; 

tip_radius = .25; 

skeleton.papillary(1).name   = 'left'; 
skeleton.papillary(1).center = left_papillary_center;
skeleton.papillary(1).radius = papillary_radius;
skeleton.papillary(1).normal = normal;

skeleton.l_to_r_papillary    = l_to_r_papillary; 

skeleton.papillary(2).name   = 'right'; 
skeleton.papillary(2).center = right_papillary_center;
skeleton.papillary(2).radius = papillary_radius;
skeleton.papillary(2).normal = normal;





