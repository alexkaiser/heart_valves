function pts = interpolate_valve_ring_points(valve, mesh, valve_ring_points_priority)

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

skeleton = valve.skeleton;

if exist('valve_ring_points_priority', 'var')
    skeleton.valve_ring_pts = valve_ring_points_priority; 
end 

if ~isfield(skeleton, 'valve_ring_pts')
    error('must have field valve.skeleton.valve_ring_pts to interpolate points'); 
end

if ~isfield(skeleton, 'ring_center')
    warning('ring_center not found, computing center');
    ring_center = mean(skeleton.valve_ring_pts,2); 
else
    ring_center = skeleton.ring_center; 
end 

if ~isfield(skeleton, 'ring_offset_angle')
    warning('setting ring_offset_angle to default value of zero'); 
    ring_offset_angle = 0; 
else 
    ring_offset_angle = skeleton.ring_offset_angle; 
end

% these are the angles at which the points lie
x = skeleton.valve_ring_pts(1,:) - ring_center(1); 
y = skeleton.valve_ring_pts(2,:) - ring_center(2); 

theta = atan2(y,x); 

[theta_unique idx] = unique(theta); 

ring_pts_unique = skeleton.valve_ring_pts(:,idx); 

mesh = mesh + ring_offset_angle; 

% center mesh on zero, because atan2 returns results on this interval 
mesh = mod(mesh, 2*pi) - pi; 

theta_three_period = [theta_unique - 2*pi, theta_unique, theta_unique + 2*pi]; 

ring_points_three_period = [ring_pts_unique, ring_pts_unique, ring_pts_unique]; 

pts = interp1(theta_three_period, ring_points_three_period', mesh, 'spline')'; 
