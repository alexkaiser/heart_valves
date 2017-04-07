function [X] = build_initial_fibers_bead_slip(leaflet)
%
% Builds initial fibers for current layout 
% 


r                       = leaflet.r; 
j_max                   = leaflet.j_max; 
k_min                   = leaflet.k_min; 
k_max                   = leaflet.k_max; 
ring_k_idx              = leaflet.ring_k_idx; 
papillary               = leaflet.papillary; 


% first and last point are in appropriate general vicinity to build initial guess 
n_papillary             = size(papillary,2); 
left_papillary          = papillary(:,1); 
right_papillary         = papillary(:,n_papillary); 

n_rings_periodic        = leaflet.n_rings_periodic; 

X = NaN * zeros(3,j_max,k_max); 

debug = false; 

if leaflet.radial_and_circumferential
    
    if n_rings_periodic > 0 
        
        % periodic, initial points on cylinder
        mesh = linspace(leaflet.min_angle, leaflet.max_angle, j_max); 
        
        for j=1:j_max
            X(:,j,ring_k_idx(j)) = [r*cos(mesh(j)); r*sin(mesh(j)); 0.0]; 
        end 
        
        % very rough physical mesh spacing
        ds = norm(X(:,1,ring_k_idx(1)) - X(:,2,ring_k_idx(2))); 
        
        for j=1:j_max 
            
            k = ring_k_idx(j) - 1;
            
            while k >= k_min(j) 
                
                % Move down on the ring 
                X(:,j,k) = X(:,j,k+1) - [0; 0; ds]; 
                k = k-1; 

            end 
            
        end         
        
    else 
    
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

