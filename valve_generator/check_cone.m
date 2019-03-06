
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

a = 1; 
r = 1.5;
h = 2; 

% limited structure here 
posterior.filter.a = a; 
posterior.filter.r = r; 
posterior.filter.h = h; 

xi  = 1.2; 
eta = 1; 

X   = cone_filter_right_cone(xi, eta, posterior)

val = cone_filter_inv_right_cone(X, posterior)

xi_calc  = val(1)
eta_calc = val(2)


N = 32; 

cone_only = false; 
if cone_only

    x = linspace(0,pi/2,N); 

    ring_quarter = [r*cos(x); r*sin(x); h*ones(size(x))]; 
    ring_plane = zeros(2,N); 

    for j=1:N
        ring_plane(:,j) = cone_filter_inv_right_cone(ring_quarter(:,j), posterior); 
    end 

    fig = figure; 
    plot(ring_plane(1,:), ring_plane(2,:)); 
    axis equal 
    title('valve ring in plane, cone only')

    mesh = linspace(0,1,N); 
    rays_in_plane = zeros(2,N,N); 

    initial_pt = [a;0]; 

    for j=1:N
        for k=1:N
            rays_in_plane(:,j,k) = mesh(j)*initial_pt  + (1-mesh(j))*ring_plane(:,k); 
        end 
    end 

    fig = figure; 
    plot(ring_plane(1,:), ring_plane(2,:));
    hold on 

    for k=1:N
        plot(rays_in_plane(1,:,k), rays_in_plane(2,:,k)); 
    end 
    axis equal 
    title('rays in plane plus valve ring, cone only')



    rays_in_3d = zeros(3,N,N); 
    for j=1:N
        for k=1:N
            rays_in_3d(:,j,k) = cone_filter(rays_in_plane(1,j,k), rays_in_plane(2,j,k), posterior); 
        end 
    end 

    fig = figure; 
    plot3(ring_quarter(1,:), ring_quarter(2,:), ring_quarter(3,:))
    hold on 
    for k=1:N
        plot3(rays_in_3d(1,:,k), rays_in_3d(2,:,k), rays_in_3d(3,:,k)); 
    end 
    axis equal
    title('cone part in space, cone only')

end 


tic 

N = 32; 
x = linspace(-3*pi/4,3*pi/4,N); 


ring_half = [r*cos(x); r*sin(x); h*ones(size(x))]; 
ring_plane = zeros(2,N); 

for j=1:N
    ring_plane(:,j) = cone_filter_inv(ring_half(:,j), posterior); 
end 

fig = figure; 
plot(ring_plane(1,:), ring_plane(2,:)); 
axis equal 
title('valve ring in plane')

mesh = linspace(0,1,N); 
rays_in_plane_right = zeros(2,N,N); 
rays_in_plane_left  = zeros(2,N,N); 

initial_pt_right = [ a;0]; 
initial_pt_left  = [-a;0]; 

for j=1:N
    for k=1:N
        rays_in_plane_right(:,j,k) = mesh(j)*initial_pt_right  + (1-mesh(j))*ring_plane(:,k); 
        rays_in_plane_left(:,j,k)  = mesh(j)*initial_pt_left   + (1-mesh(j))*ring_plane(:,k); 
    end 
end 

fig = figure; 
plot(ring_plane(1,:), ring_plane(2,:));
hold on 

for k=1:N
    plot(rays_in_plane_right(1,:,k), rays_in_plane_right(2,:,k)); 
end 
axis equal 
title('rays in plane plus valve ring')



rays_in_3d_right = zeros(3,N,N); 
rays_in_3d_left  = zeros(3,N,N); 
for j=1:N
    for k=1:N
        rays_in_3d_right(:,j,k) = cone_filter(rays_in_plane_right(1,j,k), rays_in_plane_right(2,j,k), posterior); 
        rays_in_3d_left (:,j,k) = cone_filter(rays_in_plane_left (1,j,k), rays_in_plane_left (2,j,k), posterior); 
    end 
end 

fig = figure; 
plot3(ring_half(1,:), ring_half(2,:), ring_half(3,:))
hold on 
for k=1:N
    plot3(rays_in_3d_right(1,:,k), rays_in_3d_right(2,:,k), rays_in_3d_right(3,:,k)); 
    plot3(rays_in_3d_left (1,:,k), rays_in_3d_left (2,:,k), rays_in_3d_left (3,:,k)); 
end 
axis equal
title('whole surface in space')



% invert and make sure that we get something reasonable 
rays_in_plane_right_computed = zeros(2,N,N); 
rays_in_plane_left__computed = zeros(2,N,N); 
for j=1:N
    for k=1:N
        rays_in_plane_right_computed(:,j,k) = cone_filter_inv(rays_in_3d_right(:,j,k), posterior);  
        rays_in_plane_left__computed(:,j,k) = cone_filter_inv(rays_in_3d_left (:,j,k), posterior);  
    end 
end 

fig = figure; 
plot(ring_plane(1,:), ring_plane(2,:));
hold on 

for k=1:N
    plot(rays_in_plane_right_computed(1,:,k), rays_in_plane_right_computed(2,:,k)); 
    plot(rays_in_plane_left__computed(1,:,k), rays_in_plane_left__computed(2,:,k)); 
end 
axis equal 
title('rays in plane plus valve ring - Computed by inversion')


'total time for surf build'
toc














