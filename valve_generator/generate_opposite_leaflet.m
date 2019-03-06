function leaflet = generate_opposite_leaflet(leaflet_current)
%
% Generates a leaflet which is opposite the current leaflet  
% and shares a free edge 
% 
% Free edge points are set to be NOT internal 
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

leaflet.N           = leaflet_current.N; 
leaflet.total_angle = 2*pi - leaflet_current.total_angle; 
leaflet.min_angle   = pi + leaflet.total_angle/2.0;
leaflet.max_angle   = pi - leaflet.total_angle/2.0; 

leaflet.r = leaflet_current.r; 

total_length = leaflet.total_angle * leaflet.r; 
leaflet.du = total_length / (leaflet.N+1); 
leaflet.dv = total_length / (leaflet.N+1); 

debug = false; 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
leaflet.radial_and_circumferential = leaflet_current.radial_and_circumferential; 

if ~leaflet.radial_and_circumferential
    error('generate opposite only implemented for radial and circumferential')
end 

leaflet.j_max               = leaflet_current.j_max; 
leaflet.k_max               = leaflet_current.k_max; 
leaflet.free_edge_idx_left  = leaflet_current.free_edge_idx_left;  
leaflet.free_edge_idx_right = leaflet_current.free_edge_idx_right; 
leaflet.chordae_idx_left    = leaflet_current.chordae_idx_left;  
leaflet.chordae_idx_right   = leaflet_current.chordae_idx_right; 
leaflet.is_internal         = leaflet_current.is_internal; 
leaflet.is_bc               = leaflet_current.is_bc; 
leaflet.point_idx_with_bc   = leaflet_current.point_idx_with_bc; 

% Spring constants in two directions 
leaflet.alpha    = leaflet_current.alpha; 
leaflet.beta     = leaflet_current.beta; 

% pressure takes a sign after relection 
leaflet.p_0      = -leaflet_current.p_0;


leaflet.ref_frac = leaflet_current.ref_frac; 

% No chordae here 
leaflet.chordae_tree = false; 

% Pull from other leaflet 
X = NaN * zeros(size(leaflet_current.X)); 

% Unpack, pull free edges from other leaflet 
r                       = leaflet.r; 
j_max                   = leaflet.j_max; 
k_max                   = leaflet.k_max; 
free_edge_idx_left      = leaflet.free_edge_idx_left; 
free_edge_idx_right     = leaflet.free_edge_idx_right; 

% set the valve ring
mesh = linspace(leaflet.min_angle, leaflet.max_angle, j_max); 
X(:,:,k_max) = [r*cos(mesh); r*sin(mesh); zeros(size(mesh))]; 


for i=1:size(free_edge_idx_left, 1)
    j = free_edge_idx_left(i,1); 
    k = free_edge_idx_left(i,2); 
    
    % free edges NOT internal on reflected leaflet 
    % treat as boundary condition 
    leaflet.is_internal(j,k) = false; 
    leaflet.is_bc(j,k) = true; 
end

for i=1:size(free_edge_idx_right, 1)
    j = free_edge_idx_right(i,1); 
    k = free_edge_idx_right(i,2); 
    
    % free edges NOT internal on reflected leaflet 
    % treat as boundary condition 
    leaflet.is_internal(j,k) = false; 
    leaflet.is_bc(j,k) = true; 
end 

% copy free edge 
X = sync_free_edge_to_anterior(leaflet_current, leaflet, X); 

% fix indices 
leaflet.linear_idx_offset = zeros(size(leaflet_current.linear_idx_offset)); 
count = max(max(leaflet_current.linear_idx_offset)) + 3; 
for k=1:k_max
    for j=1:j_max
        if leaflet.is_internal(j,k)
            leaflet.linear_idx_offset(j,k) = count; 
            count = count + 3; 
        end 
    end 
end


% if true, takes crude guess at closed leaflet with curvature 
pinched_interpolant = false; 

if pinched_interpolant 

    % one dimensional mesh in straight line from commissure to commissure 

    ring_l = X(:,1    ,k_max); 
    ring_r = X(:,j_max,k_max);

    line_comm_to_comm = zeros(3,k_max+1); 
    ds = 1 / (j_max - 1); 
    for m=0:j_max
        line_comm_to_comm(:,m+1) = (m*ds)*ring_r + (1 - m*ds)*ring_l; 
    end 


    % fill in fibers interpolating between free edge and ring on each side 
    for i=1:size(free_edge_idx_left, 1)

        j = free_edge_idx_left(i,1); 
        k = free_edge_idx_left(i,2); 

        % number of points on this fiber 
        num_points = k_max - k - 1; 

        % parameter spacing 
        ds = 1 / (k_max - k); 

        X_free = X(:,j,k); 
        X_ring = X(:,j,k_max); 

        X_line_commissure = 0.5 * (line_comm_to_comm(:,j) + X_free);

        num_points_first_half = floor(num_points/2); 
        ds_first = 1 / (num_points_first_half + 1); 

        num_points_second_half = num_points - num_points_first_half; 
        ds_second = 1 / (num_points_second_half + 1); 

        for m=1:num_points_first_half
            k_tmp = k + m; 
            X(:,j,k_tmp) = (m*ds_first)*X_line_commissure + (1 - m*ds_first)*X_free; 
        end

        for m=1:num_points_second_half
            k_tmp = k + m + num_points_first_half; 
            X(:,j,k_tmp) = (m*ds_second)*X_ring + (1 - m*ds_second)*X_line_commissure;                      
        end

    end

    for i=1:size(free_edge_idx_right, 1)

        j = free_edge_idx_right(i,1); 
        k = free_edge_idx_right(i,2); 

        % number of points on this fiber 
        num_points = k_max - k - 1; 

        % parameter spacing 
        ds = 1 / (k_max - k); 

        X_free = X(:,j,k); 
        X_ring = X(:,j,k_max); 

        X_line_commissure = 0.5 * (line_comm_to_comm(:,j) + X_free);

        num_points_first_half = floor(num_points/2); 
        ds_first = 1 / (num_points_first_half + 1); 

        num_points_second_half = num_points - num_points_first_half; 
        ds_second = 1 / (num_points_second_half + 1); 

        for m=1:num_points_first_half
            k_tmp = k + m; 
            X(:,j,k_tmp) = (m*ds_first)*X_line_commissure + (1 - m*ds_first)*X_free; 
        end

        for m=1:num_points_second_half
            k_tmp = k + m + num_points_first_half; 
            X(:,j,k_tmp) = (m*ds_second)*X_ring + (1 - m*ds_second)*X_line_commissure;                      
        end

    end  


    % linear interpolant 
    else 

        for i=1:size(free_edge_idx_left, 1)

            j = free_edge_idx_left(i,1); 
            k = free_edge_idx_left(i,2); 

            % number of points on this fiber 
            num_points = k_max - k - 1; 

            % parameter spacing 
            ds = 1 / (k_max - k); 

            X_free = X(:,j,k); 
            X_ring = X(:,j,k_max); 

            for m=1:num_points
                k_tmp = k + m; 
                X(:,j,k_tmp) = (m*ds)*X_ring + (1 - m*ds)*X_free; 
            end 

        end


        for i=1:size(free_edge_idx_right, 1)

            j = free_edge_idx_right(i,1); 
            k = free_edge_idx_right(i,2); 

            % number of points on this fiber 
            num_points = k_max - k - 1; 

            % parameter spacing 
            ds = 1 / (k_max - k); 

            X_free = X(:,j,k); 
            X_ring = X(:,j,k_max); 

            for m=1:num_points
                k_tmp = k + m; 
                X(:,j,k_tmp) = (m*ds)*X_ring + (1 - m*ds)*X_free; 
            end

        end 

end 

leaflet.X = X; 


% NaN mask so using bad values will give errors 
leaflet.R = NaN * zeros(size(leaflet.X)); 

for i=1:size(free_edge_idx_left, 1)
    j = free_edge_idx_left(i,1); 
    k = free_edge_idx_left(i,2); 
    
    % Free edge and neighbors have rest positions 
    % Set to current position for now 
    % leaflet.R(:,j  ,k  ) = leaflet.X(:,j  ,k  );  
    leaflet.R(:,j+1,k  ) = leaflet.X(:,j+1,k  );  
    leaflet.R(:,j  ,k+1) = leaflet.X(:,j  ,k+1);  
    
    leaflet.X(:,j  ,k  ) = NaN * leaflet.X(:,j  ,k  ); 
end

for i=1:size(free_edge_idx_right, 1)
    j = free_edge_idx_right(i,1); 
    k = free_edge_idx_right(i,2); 
   
    % right free edge has neighbors up in k 
    % but down in j
    % leaflet.R(:,j  ,k  ) = leaflet.X(:,j  ,k  );  
    leaflet.R(:,j-1,k  ) = leaflet.X(:,j-1,k  );  
    leaflet.R(:,j  ,k+1) = leaflet.X(:,j  ,k+1);
    
    leaflet.X(:,j  ,k  ) = NaN * leaflet.X(:,j  ,k  ); 
end 


if debug 
    figure; 

    x_component = squeeze(X(1,:,:)); 
    y_component = squeeze(X(2,:,:)); 
    z_component = squeeze(X(3,:,:)); 

    width = 1.5; 
    surf(x_component, y_component, z_component, 'LineWidth',width);

    axis equal 
    axis auto 

    xlabel('x'); 
    ylabel('y'); 
    title('leaflet')
end 




