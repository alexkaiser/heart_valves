function [X, commissure_points, center, comm_position_set] = interpolate_valve_ring_pts_aortic(leaflet, valve_ring_pts, height_min_comm)


j_max = leaflet.j_max; 
k_max = leaflet.k_max; 

X = NaN * zeros(3,j_max,k_max); 



[m n_ring_pts] = size(valve_ring_pts); 

if (m ~= 3)
    error('valve_ring_points must be 3 by n_ring_pts'); 
end 

interp_points_initial = linspace(0,1,n_ring_pts);

upsample_multiplier = 1000; 
n_pts_upsampled = n_ring_pts*upsample_multiplier;
interp_points_upsampled = linspace(0,1,n_pts_upsampled);

valve_ring_upsampled = interp1(interp_points_initial, valve_ring_pts', interp_points_upsampled, 'pchip')'; 


debug = false; 
if debug 
    figure; 
    plot3(valve_ring_upsampled(1,:), valve_ring_upsampled(2,:), valve_ring_upsampled(3,:))
    hold on 
    plot3(valve_ring_pts(1,:), valve_ring_pts(2,:), valve_ring_pts(3,:), '*')

end 

% compute cumumlative differences 
differences = valve_ring_upsampled(:,2:end) - valve_ring_upsampled(:,1:(end-1));
norm_diffs = vecnorm(differences);

cumulative_distance = cumsum(norm_diffs);
total_arc_length = cumulative_distance(end);

l_comm_min_range = 1; 
l_comm_max_range = find(cumulative_distance > height_min_comm, 1, 'first');

r_comm_min_range = find(cumulative_distance > (total_arc_length - height_min_comm), 1, 'first'); 
r_romm_max_range = n_pts_upsampled;


X_left_comm_interp = valve_ring_upsampled(:, l_comm_min_range:l_comm_max_range); 
% these go from bottom to top 
X_left_comm_interp = flip(X_left_comm_interp,2);

X_annulus_bottom_interp = valve_ring_upsampled(:, l_comm_max_range:r_comm_min_range); 
X_right_comm_interp = valve_ring_upsampled(:, r_comm_min_range:r_romm_max_range); 

% left comm points 
X(:,1,:) = compute_equally_spaced_pts(X_left_comm_interp, k_max);
X_l_comm_min = X(:,1,1);

% right comm points 
X(:,j_max,:) = compute_equally_spaced_pts(X_right_comm_interp, k_max);
X_r_comm_min = X(:,j_max,1); 

% bottom annulus points 
% overwrite corner points (which should be equal)
X(:,:,1) = compute_equally_spaced_pts(X_annulus_bottom_interp, j_max);


% error checking on positions
tol = 1e-14;
if norm(X_l_comm_min - X(:,1,1)) > tol 
    warning('inconsistent position at base of radial attachment')
end 
if norm(X_r_comm_min - X(:,j_max,1)) > tol
    warning('inconsistent position at base of radial attachment')
end 

% commissure points
% just the first and last point on the free edge 
commissure_points = [X(:,1,k_max), X(:,j_max,k_max)];

% vector mean is middle interpolation point 
center = mean(commissure_points, 2);

comm_position_set = true; 
  
debug_spacing = false; 
if debug_spacing
    figure; 

    x_component = squeeze(X(1,:,1)); 
    y_component = squeeze(X(2,:,1)); 
    z_component = squeeze(X(3,:,1)); 

    hold on 
    width = 1.0; 
    plot3(x_component, y_component, z_component, 'ko', 'LineWidth',width);

    x_component_valve_ring_upsampled = squeeze(valve_ring_upsampled(1,:)); 
    y_component_valve_ring_upsampled = squeeze(valve_ring_upsampled(2,:)); 
    z_component_valve_ring_upsampled = squeeze(valve_ring_upsampled(3,:)); 
    
    plot3(x_component_valve_ring_upsampled, y_component_valve_ring_upsampled, z_component_valve_ring_upsampled, 'r', 'LineWidth',width);
    
    x_component = squeeze(X(1,1,:)); 
    y_component = squeeze(X(2,1,:)); 
    z_component = squeeze(X(3,1,:)); 
    plot3(x_component, y_component, z_component, '*', 'LineWidth',width);
    
    x_component = squeeze(X(1,j_max,:)); 
    y_component = squeeze(X(2,j_max,:)); 
    z_component = squeeze(X(3,j_max,:)); 
    plot3(x_component, y_component, z_component, '+', 'LineWidth',width);

    axis equal 
    axis auto 
    
    xlabel('x'); 
    ylabel('y'); 
    title('spacing ')
    
    % quick length consistency check 
    tol_lengths = 1e-3; 
    max_diff_absolute = 0; 
    
    
    % comm point below 
    arc_len_first = norm(X(:,1,1) - X(:,2,1)); 

    for j=2:j_max
        j
        arc_len_current = norm(X(:,j-1,1) - X(:,j,1)) 
        max_diff_absolute = max(max_diff_absolute, abs((arc_len_current - arc_len_first)));
        rel_err = abs((arc_len_current - arc_len_first)/arc_len_first)
        if rel_err > tol_lengths
            warning('inconsistent arc lengths after equalizing');                     
        end 
    end 

    max_diff_absolute

end 





