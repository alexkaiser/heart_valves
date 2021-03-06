function [X] = build_initial_fibers_bead_slip(leaflet, valve)
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


r                = leaflet.r; 
j_max            = leaflet.j_max; 
k_min            = leaflet.k_min; 
k_max            = leaflet.k_max; 
ring_k_idx       = leaflet.ring_k_idx; 
n_rings_periodic = leaflet.n_rings_periodic; 

N_anterior       = valve.N_anterior; 
N_posterior      = valve.N_posterior;

commissural_leaflets = valve.commissural_leaflets; 
if commissural_leaflets
    N_commissure = valve.N_commissure; 
end 
    
X = NaN * zeros(3,j_max,k_max); 

debug = false; 

if leaflet.radial_and_circumferential
    
    if true % n_rings_periodic > 0 
        
        if commissural_leaflets
            
            % centered around zero 
            min_anterior = -leaflet.total_angle_anterior/2; 
            max_anterior =  leaflet.total_angle_anterior/2; 

            % mesh anterior inclusive of ends 
            mesh_anterior   = linspace(min_anterior, max_anterior, N_anterior); 

            min_posterior = pi - leaflet.total_angle_posterior/2; 
            max_posterior = pi + leaflet.total_angle_posterior/2;             
            
            % calculate inclusive, the crop outer points 
            mesh_right_comm = linspace(max_anterior, min_posterior, N_commissure + 2); 
            mesh_right_comm = mesh_right_comm(2:(end-1)); 
            
            mesh_posterior  = linspace(min_posterior, max_posterior, N_posterior);  
            
            min_anterior_wrapped = min_anterior + 2*pi; 
            mesh_left_comm  = linspace(max_posterior, min_anterior_wrapped, N_commissure + 2);
            mesh_left_comm  = mesh_left_comm(2:(end-1)); 
            
            mesh = [mesh_anterior mesh_right_comm mesh_posterior mesh_left_comm];
            
        else 
            
            % centered around zero 
            min_anterior = -leaflet.total_angle_anterior/2; 
            max_anterior =  leaflet.total_angle_anterior/2; 

            % mesh anterior inclusive of ends 
            mesh_anterior   = linspace(min_anterior, max_anterior, N_anterior); 

            % mesh posterior includes two anterior points
            min_anterior_wrapped = min_anterior + 2*pi; 
            mesh_posterior  = linspace(max_anterior, min_anterior_wrapped, N_posterior + 2);
            mesh_posterior  = mesh_posterior(2:(end-1)); 

            mesh = [mesh_anterior mesh_posterior];
        end 
        
        if isfield(valve.skeleton, 'valve_ring_pts')
            % use whatever points are measured for valve ring 
            for j=1:j_max
                X(:,j,ring_k_idx(j)) = interpolate_valve_ring_points(valve, mesh(j)); 
                % [r*(cos(mesh(j)) + x_coord_extra(mesh(j))); r*sin(mesh(j)); 0.0]; 
            end 

        else
            
            % simple analytic shape for valve ring 
            if isfield(valve, 'dip_anterior_systole') && valve.dip_anterior_systole 
                x_coord_extra = @(t) cos_bump(t, valve.total_angle_dip, valve.r_dip); 
            else 
                x_coord_extra = @(t) 0 .* t; 
            end     


            for j=1:j_max
                X(:,j,ring_k_idx(j)) = [r*(cos(mesh(j)) + x_coord_extra(mesh(j))); r*sin(mesh(j)); 0.0]; 
            end 
        
        end
        
        % very rough physical mesh spacing
        % should be unimportant since this is just for an initial guess 
        ds = 0.5 * norm(X(:,1,ring_k_idx(1)) - X(:,2,ring_k_idx(2))); 
        
        for j=1:j_max 
            
            k = ring_k_idx(j) - 1;
            
            while k >= k_min(j) 
                
                % Move down on the ring 
                X(:,j,k) = X(:,j,k+1) - [0; 0; ds]; 
                k = k-1; 

            end 
            
        end         
        
    else 
        
        % first and last point are in appropriate general vicinity to build initial guess 
        papillary               = leaflet.papillary; 
        n_papillary             = size(papillary,2);    
        left_papillary          = papillary(:,1);
        right_papillary         = papillary(:,n_papillary); 
        
        % set the valve ring

        % previous version
        % this keeps the endpoints fixed, but the internal points do not line up 
        % when the mesh is refined 
        mesh = linspace(leaflet.min_angle, leaflet.max_angle, j_max); 

        % one extra point, 
    %     mesh = linspace(leaflet.min_angle, leaflet.max_angle, j_max + 1);
    %     
    %     % which is then thrown out
    %     % take the rightmost point on the anterior 
    %     % leftmost on posterior, arbitrarily 
    %     if leaflet.min_angle < leaflet.max_angle
    %         mesh = mesh(1:j_max); 
    %     else 
    %         mesh = mesh(2:(j_max+1)); 
    %     end 

        % X(:,:,k_max) = [r*cos(mesh); r*sin(mesh); zeros(size(mesh))]; 

        for j=1:j_max
            X(:,j,ring_k_idx(j)) = [r*cos(mesh(j)); r*sin(mesh(j)); 0.0]; 
        end 

        % Set free edge according to interpolating surface

        ring_l = X(:, 1    , ring_k_idx(1    )); 
        ring_r = X(:, j_max, ring_k_idx(j_max)); 

        ds = 1/(j_max - 1); 

        interpolating_surf = @(s,t) t*(s*ring_r + (1-s)*ring_l) + (1-t)*(s*right_papillary + (1-s)*left_papillary);   
        t_of_s = @(s) abs(s-1/2) + 1/2 - ds; 
    
        for j=1:j_max 
            k = k_min(j); 

            s = (j-1)*ds; 
            X(:,j,k) = interpolating_surf(s,t_of_s(s));
        end     

        if debug 
            figure; 
            plot3(X(1,:), X(2,:), X(3,:), '-o'); 
            xlabel('x'); 
            ylabel('y'); 
            title('free edge and reference')
        end 

        % linear interpolant from free edge 
        % fill in fibers interpolating between free edge and ring 
        for j=1:j_max 
            k = k_min(j); 

            % number of points on this fiber 
            num_points = ring_k_idx(j) - k - 1; 

            % parameter spacing 
            ds = 1 / (ring_k_idx(j) - k); 

            X_free = X(:,j,k); 
            X_ring = X(:,j,ring_k_idx(j)); 

            for m=1:num_points
                k_tmp = k + m; 
                X(:,j,k_tmp) = (m*ds)*X_ring + (1 - m*ds)*X_free; 
            end 

        end 
    
    end 
    
%     % translate solution if center of the ring is not the origin
%     % as it would be in the general case 
%     if isfield(valve.skeleton, 'ring_center')
%         X = X + valve.skeleton.ring_center; 
%     end 
    
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
    

else 
    error('diag not implemented with bead slip'); 
end 

