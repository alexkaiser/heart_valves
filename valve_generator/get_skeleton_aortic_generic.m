function skeleton = get_skeleton_aortic_generic(r, h1, hc, r_commissure)
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

% https://e-echocardiography.com/page/page.php?UID=1867001
% middle of range, simple hack 
% skeleton.r = 1.0; 
% skeleton.r = 1.25; 
% skeleton.r_commissure = 1.0 * skeleton.r; 
% 
% skeleton.normal_height = 1.4 * skeleton.r; 
% 
% skeleton.height_min_comm = 0.6 * skeleton.r; 
% 
% skeleton.ring_offset_angle = 0; 

% from... 
% A general three-dimensional parametric geometry of the native aortic valve and root for biomechanical modeling
% 
% approx radius from normal 1 
% .05 larger than value in paper 
if ~exist('r','var')
    r = 1.25; 
end 

if ~exist('r_commissure','var')
    r_commissure = r; 
end 

skeleton.r = r; 

% aorta radius 
% r_aorta = r; % 1.1 * r; 

% sinus radius at widest point (valve not to here)
% r_sinus = r; %1.4 * r; 

% commissure radius 
skeleton.r_commissure = r_commissure; 

% r_co in paper 
% wide point in verticle plane through origin and commissure 
% skeleton.r_co = (1/2) * (r_aorta + r_sinus); 

% height from annulus to origin 
% equal to height from annulus to height at which r_co is the radius
if ~exist('h1','var')
    h1 = .9 * r; 
end 

% commissure height measured from origin (annulus is not at the origin)
if ~exist('hc','var')
    hc = .5 * r; 
end 

skeleton.height_min_comm = h1; 

% totall commissures 
skeleton.normal_height = h1 + hc; 

skeleton.ring_offset_angle = 0; 



heights = [0; skeleton.height_min_comm; skeleton.normal_height]; 
radii   = [r; skeleton.r_commissure; skeleton.r_commissure]; 
skeleton.r_of_z = @(z) (z<0) .* r + ... 
                       (0<=z).*(z<skeleton.normal_height) .* interp1(heights, radii, z, 'pchip') + ... 
                       (skeleton.normal_height<=z) .* skeleton.r_commissure; 


                                    
if debug 
    r_of_z_spline = @(z) (z<0) .* r + ... 
                       (0<=z).*(z<skeleton.normal_height) .* interp1(heights, radii, z, 'spline') + ... 
                       (skeleton.normal_height<=z) .* skeleton.r_commissure; 
    
    r_of_z_linear = @(z) (z<0) .* r + ... 
       (0<=z).*(z<skeleton.normal_height) .* interp1(heights, radii, z, 'linear') + ... 
       (skeleton.normal_height<=z) .* skeleton.r_commissure; 
    
    z = linspace(0 - 0.02, skeleton.normal_height + .02, 1000); 
    plot(skeleton.r_of_z(z), z); 
    hold on 
    plot(r_of_z_spline(z), z); 
    plot(r_of_z_linear(z), z); 
    legend('pchip', 'spline', 'linear')
    
    xlabel('r_of_z')
    ylabel('z')
    title('spline for radius')
%     xlim([0 1.5])
%     ylim([0 skeleton.normal_height])
    axis equal
        
end 


