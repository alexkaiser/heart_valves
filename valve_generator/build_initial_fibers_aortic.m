function [X] = build_initial_fibers_aortic(leaflet, valve)
%
% Builds initial fibers for current layout 
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

j_max            = leaflet.j_max; 
k_max            = leaflet.k_max; 

N = leaflet.N; 
N_each = leaflet.N_each; 
    
X = NaN * zeros(3,j_max,k_max); 

debug = false; 

free_edge_smooth = false; 

if isfield(valve.skeleton, 'valve_ring_pts')
    % use whatever points are measured for valve ring 
%     for j=1:j_max
%         X(:,j,ring_k_idx(j)) = interpolate_valve_ring_points(valve, mesh(j)); 
%         % [r*(cos(mesh(j)) + x_coord_extra(mesh(j))); r*sin(mesh(j)); 0.0]; 
%     end 

    error('arbitrary skeleton not implemented for aortic')
else

    % circular rings plus interpolants 
    
    r                 = valve.skeleton.r; 
    r_commissure      = valve.skeleton.r_commissure; 
    normal_height     = valve.skeleton.normal_height; 
    ring_offset_angle = valve.skeleton.ring_offset_angle; 
    
    % valve ring lives at k=1
    
    % cusp radius is one sixth of the circumference 
    free_edge_cusp_radius = 2*pi*r/6; 
    
    if free_edge_cusp_radius > normal_height
        error('inconsistent radius and height'); 
    end 
    
    du = leaflet.du; 
    if abs(du - 1/N) > eps
        error('inconsistent values in mesh'); 
    end 
    
    % minimum ring 
    k = 1; 
    for j = 1:j_max
        
        j_this_cusp = mod(j,N_each); 
        x_this_cusp = j_this_cusp / N_each; 
        center_cusp = 1/2; 
        
        circle_height = free_edge_cusp_radius * 2 * sqrt((1/2)^2 - (x_this_cusp - center_cusp)^2); 
        
        % top of circular part (bottom of commissure) minus a circle 
        z_tmp = free_edge_cusp_radius - circle_height; 
        
        X(:,j,k) = r * [cos(j*du*2*pi + ring_offset_angle) ; sin(j*du*2*pi + ring_offset_angle); z_tmp]; 
    end
    
    % commissure points
    commissure_points = zeros(3,3); 
    for j = (1:3)
        commissure_points(:,j) = [r_commissure * cos((j-1)*N_each*du*2*pi + ring_offset_angle); ...
                                  r_commissure * sin((j-1)*N_each*du*2*pi + ring_offset_angle); ...
                                  normal_height]; 
    end
     
    % center of commissure points 
    center = [0; 0; free_edge_cusp_radius + 0.5 * (normal_height - free_edge_cusp_radius)]; 
    
    k = k_max; 
    
    for comm_idx = 1:3
        
        % point one internal of commissure to point that m
        % N_each is a power of two 
        min_idx = (comm_idx-1)*N_each;         
        
        dj_interp = 1/(N_each/2); 
        
        for j=1:(N_each/2)
            X(:,j + min_idx           ,k) = (1 - j*dj_interp) * commissure_points(:,comm_idx) + j*dj_interp * center; 
        end 
        
        % center out to commissure point 
        if comm_idx < 3
            comm_idx_next = comm_idx+1; 
        else 
            comm_idx_next = 1; 
        end
        
        for j=1:(N_each/2)
            X(:,j + min_idx + N_each/2,k) = (1 - j*dj_interp) * center + j*dj_interp * commissure_points(:,comm_idx_next); 
        end 
        
        if free_edge_smooth && (N_each > 16)
            
            smooth_points = N_each/16; 
            
            free_edge_copy = X(:,:,k); 
            
            smooth_min = (N_each/2) - smooth_points; 
            smooth_max = (N_each/2) + smooth_points; 
            
            for j = smooth_min:smooth_max
                for component = 1:3
                    free_edge_copy(component,j + min_idx) = mean(X(component,(j + min_idx - smooth_points):(j + min_idx + smooth_points)  ,k)); 
                end 
            end 
            
            X(:,:,k) = free_edge_copy; 
            
        end 

    end 
        
%     % now, linear interpolation from annulus to free edge 
    for j=1:j_max 
        dk_interp = 1/(k_max-1);         
        for k=2:(k_max-1)
            X(:,j,k) = (1 - (k-1)*dk_interp) * X(:,j,1) + (k-1)*dk_interp * X(:,j,k_max); 
        end         
    end 

end



if debug 
    figure; 

    x_component = squeeze(X(1,:,:)); 
    y_component = squeeze(X(2,:,:)); 
    z_component = squeeze(X(3,:,:)); 

    width = 1.0; 
    surf(x_component, y_component, z_component, 'LineWidth',width);

%     hold on 
%     plot3(x_component, y_component, z_component, 'ko', 'LineWidth',width);
    
    axis equal 
    axis auto 

    xlabel('x'); 
    ylabel('y'); 
    title('aortic leaflets initial')
    
    figure; 
    % six times circle radius
    z_free_edge = squeeze(X(3,:,1));
    x_circles = linspace(0,6 * free_edge_cusp_radius, length(z_component)); 
    plot(x_circles, z_free_edge, 'ko'); 
    axis equal 
    title('circles in 2d')
    
end 
    

