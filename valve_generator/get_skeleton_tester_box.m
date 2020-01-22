function skeleton = get_skeleton_tester_box()
% 
% hardcoded patient specific MV skeleton
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

debug = false; 

inches_to_cm = 2.54; 

% ring_points_x = [-0.561
% -0.541
% -0.405
% -0.184
% 0.116
% 0.361
% 0.498
% 0.576
% 0.588
% -0.541
% -0.405
% -0.184
% 0.116
% 0.361
% 0.498
% 0.576]; 
% 
% ring_points_y = [0
% -0.199
% -0.499
% -0.699
% -0.787
% -0.699
% -0.499
% -0.199
% 0
% 0.199
% 0.499
% 0.699
% 0.787
% 0.699
% 0.499
% 0.199]; 
% 
% ring_points_z = [0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27
% 0.27]; 

% ring_points_inches = [ring_points_x'; ring_points_y'; ring_points_z']; 

ring_points_inches = [
    0	-0.561	0.27
-0.199	-0.541	0.27
-0.499	-0.405	0.27
-0.699	-0.184	0.27
-0.787	0.116	0.27
-0.699	0.361	0.27
-0.499	0.498	0.27
-0.199	0.576	0.27
0	0.588	0.27
0.199	-0.541	0.27
0.499	-0.405	0.27
0.699	-0.184	0.27
0.787	0.116	0.27
0.699	0.361	0.27
0.499	0.498	0.27
0.199	0.576	0.27]'; 

ring_z = ring_points_inches(3,1); 

ring_points_inches(3,:) = ring_points_inches(3,:) - ring_z;

if any(ring_points_inches(3,:) > eps)
    error('should have zero z coordinate'); 
end 

valve_ring_pts = inches_to_cm * ring_points_inches; 

skeleton.valve_ring_pts = valve_ring_pts; 

skeleton.n_ring_pts = size(skeleton.valve_ring_pts,2); 

% origin 
skeleton.ring_center = [0;0;0]; %mean(skeleton.valve_ring_pts')'; 

radii = zeros(skeleton.n_ring_pts, 1); 
for j=1:skeleton.n_ring_pts
    radii(j) = norm(skeleton.valve_ring_pts(:,j) - skeleton.ring_center); 
end

skeleton.r = mean(radii); 


tip_radius = .2; 

papillary_left   = [.787, .116, -1.23] + [0 0 -ring_z];
papillary_left   = papillary_left * inches_to_cm; 
papillary_left   = papillary_left'; 
papillary_right  = [-.787, .116, -1.23] + [0 0 -ring_z];
papillary_right  = papillary_right * inches_to_cm;  
papillary_right  = papillary_right';



skeleton.ring_offset_angle = 3*pi/2; 

l_to_r_papillary = papillary_right - papillary_left; 
l_to_r_papillary = l_to_r_papillary / norm(l_to_r_papillary);

skeleton.papillary(1).name   = 'left'; 
skeleton.papillary(1).center = papillary_left;
skeleton.papillary(1).radius = tip_radius;
skeleton.papillary(1).normal = [0; 0; 1];

skeleton.l_to_r_papillary    = l_to_r_papillary; 

skeleton.papillary(2).name   = 'right'; 
skeleton.papillary(2).center = papillary_right;
skeleton.papillary(2).radius = tip_radius;
skeleton.papillary(2).normal = [0; 0; 1];


if debug 
    fig = figure; 
    
    x = squeeze(skeleton.valve_ring_pts(1,:));
    y = squeeze(skeleton.valve_ring_pts(2,:));
    z = squeeze(skeleton.valve_ring_pts(3,:));
    plot3(x,y,z,'k'); 
    hold on 
    plot3(papillary_left(1), papillary_left(2), papillary_left(3), 'ko')
    plot3(papillary_right(1), papillary_right(2), papillary_right(3), 'k*')
    axis equal 
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('valve ring tester box')
end 





