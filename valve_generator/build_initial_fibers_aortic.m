function [X, len_annulus_min, len_annulus_each] = build_initial_fibers_aortic(leaflet, valve, height_min_comm_override, power_override)
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

j_max = leaflet.j_max; 
k_max = leaflet.k_max; 

N = leaflet.N; 
N_each = leaflet.N_each; 

if isfield(leaflet, 'N_leaflets')
    N_leaflets = leaflet.N_leaflets; 
else 
    % default tricuspid
    N_leaflets = 3; 
end 
X = NaN * zeros(3,j_max,k_max); 

debug = true; 


if isfield(valve.skeleton, 'valve_ring_pts')
    % use whatever points are measured for valve ring 
 
    [pts, commissure_points, center] = interpolate_valve_ring_pts_aortic(valve, N_each, N_leaflets); 
    X(:,:,1) = pts; 
    % error('arbitrary skeleton not implemented for aortic')
else

    % circular rings plus interpolants 
    
    r                 = valve.skeleton.r; 
    r_commissure      = valve.skeleton.r_commissure; 
    normal_height     = valve.skeleton.normal_height;
    ring_offset_angle = valve.skeleton.ring_offset_angle; 
    
    % valve ring lives at k=1
    
    % cusp radius is one sixth of the circumference 
    free_edge_cusp_radius = 2*pi*r/6; 
    
    if isfield(valve.skeleton, 'height_min_comm')
        height_min_comm = valve.skeleton.height_min_comm; 
    else 
        height_min_comm = free_edge_cusp_radius; 
    end 
    
    if exist('height_min_comm_override', 'var')
        % override existing height_min_comm
        
        height_comm = normal_height - height_min_comm; 
        
        height_min_comm = height_min_comm_override; 
        normal_height = height_min_comm + height_comm; 
    end 
    
    % polynomial height profile 
    % default value 
    power = 3; 
    if isfield(valve, 'annulus_power')
        power = valve.annulus_power; 
    end 
    if exist('power_override', 'var')
        power = power_override; 
        if isfield('valve', 'annulus_power')
            warning('annulus_power and power_override both provided, using power_override')
        end 
    end 

    if isfield(valve, 'use_annulus_flattened_pts') && valve.use_annulus_flattened_pts 
        
        if ~isfield(valve, 'annulus_flattened_normalized')
            error('annulus_flattened_normalized required if valve.use_annulus_flattened_pts is true'); 
        end 
        
        [m,n] = size(valve.annulus_flattened_normalized); 
        if n ~= 2
            error('valve.annulus_flattened_normalized must be two column matrix')
        end 
        
        if (abs(valve.annulus_flattened_normalized(1,2) - 1) > eps) || ... 
           (abs(valve.annulus_flattened_normalized(m,2) - 1) > eps) || ... 
           (abs(valve.annulus_flattened_normalized(1,1) ~= 0)) || ... 
           (abs(valve.annulus_flattened_normalized(m,1) ~= 1)) 
            error('valve.annulus_flattened_normalized boundary values incorrect')
        end 
        
        z_tmp_fn = @(x_this_cusp) height_min_comm * interp1(valve.annulus_flattened_normalized(:,1), valve.annulus_flattened_normalized(:,2), x_this_cusp, 'pchip'); 
    else
        % default polynomial 
        center_cusp = 1/2;
        normalization = 2^power; % normalization so function takes value 1 at 1/2 
        z_tmp_fn = @(x_this_cusp) height_min_comm * normalization * abs(x_this_cusp - center_cusp)^power; 
    end 
    
    du = leaflet.du; 
    if abs(du - 1/N) > eps
        error('inconsistent values in mesh'); 
    end 
    
    % annulus_points_even_spacing = true; 
    
    if isfield(valve, 'annulus_points_even_spacing') && valve.annulus_points_even_spacing
        error('not implemented')
        % just put a bunch of points 
        % then interpolate to equally spaced 
        mesh_scaling = 1000; 
        N_interp = N * mesh_scaling; 
        N_each_interp = N_each * mesh_scaling; 
        du_interp = 1/N_interp; 
        j_max_interp = N_each_interp + 1; 
        X_annulus_interp = zeros(3,j_max_interp); 
        
        % initial positions to interpolate from 
        for j = 1:(N_each_interp + 1)

            x_this_cusp = (j-1) / N_each_interp; 

            z_tmp = z_tmp_fn(x_this_cusp);  

            r_tmp = r; 

            if r ~= r_commissure 
                if isfield(valve.skeleton, 'r_of_z')
                    r_tmp = valve.skeleton.r_of_z(z_tmp); 
                else 
                    error('this needs to compute r(z) here if r ~= r_commissure'); 
                end 
            end         

            X_annulus_interp(:,j) = [r_tmp*cos((j-1)*du_interp*2*pi + ring_offset_angle) ; r_tmp*sin((j-1)*du_interp*2*pi + ring_offset_angle); z_tmp]; 
        end        
        
        arc_len_cumulative = zeros(N_each_interp + 1,1); 

        for j=2:(N_each_interp+1)
            arc_len_cumulative(j) = arc_len_cumulative(j-1) + norm(X_annulus_interp(:,j-1) - X_annulus_interp(:,j)); 
        end 

        arc_len_total = arc_len_cumulative(end); 

        query_pts = (arc_len_total/N_each) * (0:N_each);  

        X(:,1:j_max,1)= interp1(arc_len_cumulative, X_annulus_interp', query_pts)'; 
            
                         
        
        debug_spacing = true; 
        if debug_spacing
            figure; 

            x_component = squeeze(X(1,:,1)); 
            y_component = squeeze(X(2,:,1)); 
            z_component = squeeze(X(3,:,1)); 

            hold on 
            width = 1.0; 
            plot3(x_component, y_component, z_component, 'ko', 'LineWidth',width);

            x_component_X_annulus_interp = squeeze(X_annulus_interp(1,:)); 
            y_component_X_annulus_interp = squeeze(X_annulus_interp(2,:)); 
            z_component_X_annulus_interp = squeeze(X_annulus_interp(3,:)); 
            
            plot3(x_component_X_annulus_interp, y_component_X_annulus_interp, z_component_X_annulus_interp, 'r', 'LineWidth',width);
            
            axis equal 
            axis auto 
            
            xlabel('x'); 
            ylabel('y'); 
            title('spacing ')
            
            % quick length consistency check 
            tol_lengths = 1e-3; 
            
            % comm point below 
            arc_len_first = norm(X(:,1,1) - X(:,2,1)); 

            for j=2:N_each
                arc_len_current = norm(X(:,j-1,1) - X(:,j,1)); 

                if abs((arc_len_current - arc_len_first)/arc_len_first) > tol_lengths
                    error('inconsistent arc lengths after equalizing'); 
                end 
            end 
                                              
        end 
                
        
        
    else 
        % minimum ring 
        k = 1; 
        for j = 1:j_max
            
            % no mod necessary 
            % j_this_cusp = mod(j,N_each); 
            % subtract one for bc
            x_this_cusp = (j-1) / N_each; 

            z_tmp = z_tmp_fn(x_this_cusp); 

            r_tmp = r; 

            if r ~= r_commissure 
                if isfield(valve.skeleton, 'r_of_z')
                    r_tmp = valve.skeleton.r_of_z(z_tmp); 
                else 
                    error('this needs to compute r(z) here if r ~= r_commissure'); 
                end 
            end         


            X(:,j,k) = [r_tmp*cos((j-1)*du*2*pi + ring_offset_angle) ; r_tmp*sin((j-1)*du*2*pi + ring_offset_angle); z_tmp]; 
        end
    end 
    
    if debug 
        annulus_z = X(3,:,1); 
        max_z_annulus = max(annulus_z); 
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

end 


% set initial reasonable free edge 
% for interpolated or parameteric 
k = k_max; 
comm_idx = 1;
dj_interp = 1/(N_each/2); 

% first half of leaflet including commissure 
for j=0:(N_each/2)
    X(:,j + 1,k) = (1 - j*dj_interp) * commissure_points(:,comm_idx) + j*dj_interp * center; 
end 

% center out to commissure point 
if comm_idx < N_leaflets
    comm_idx_next = comm_idx+1; 
else 
    comm_idx_next = 1; 
end

% second half of leaflet 
for j=1:(N_each/2)
    X(:,j + 1 + N_each/2,k) = (1 - j*dj_interp) * center + j*dj_interp * commissure_points(:,comm_idx_next); 
end 



% linear interpolation from annulus to free edge 
for j=1:j_max 
    dk_interp = 1/(k_max-1);         
    for k=2:(k_max-1)
        X(:,j,k) = (1 - (k-1)*dk_interp) * X(:,j,1) + (k-1)*dk_interp * X(:,j,k_max);  
    end         
end 


% minimum ring 
len_annulus_min = 0; 
len_annulus_each = zeros(j_max,1); 
k = 1; 
for j = 1:N_each

    [valid, j_nbr, k_nbr] = get_indices(leaflet, j, k, j+1, k); 
    if ~valid
        error('failed to find valid index'); 
    end 

    len_annulus_min = len_annulus_min + norm(X(:,j_nbr,k_nbr) - X(:,j,k));     
    len_annulus_each(j) = norm(X(:,j_nbr,k_nbr) - X(:,j,k)); 
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
    
%     figure; 
%     % six times circle radius
%     z_free_edge = squeeze(X(3,:,1));
%     
%     dx = 6 * free_edge_cusp_radius / length(z_component); 
%     x_circles = dx * (1:length(z_component));     
%     % x_circles = linspace(0,6 * free_edge_cusp_radius, length(z_component)); 
%     plot(x_circles, z_free_edge, 'ko'); 
%     axis equal 
%     title('circles in 2d')
%     hold on 
%     
%     th = 0:.001:2*pi;         
%     x = free_edge_cusp_radius * cos(th) + free_edge_cusp_radius; 
%     y = free_edge_cusp_radius * sin(th) + free_edge_cusp_radius; 
%     plot(x,y); 
    
end 
    

