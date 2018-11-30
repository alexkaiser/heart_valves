function pts = interpolate_valve_ring_points(valve, mesh)

skeleton = valve.skeleton;

if ~isfield(skeleton, 'valve_ring_pts')
    error('must have field valve.skeleton.valve_ring_pts to interpolate points'); 
end

if ~isfield(skeleton, 'ring_center')
    warning('ring_center not found, using origin');
    ring_center = zeros(3,1); 
else
    ring_center = skeleton.ring_center; 
end 

if ~isfield(skeleton, 'ring_offet_angle')
    warning('setting ring_offset_angle to default value of zero'); 
    ring_offet_angle = 0; 
else 
    ring_offet_angle = skeleton.ring_offset_angle; 
end

% these are the angles at which the points lie
x = skeleton.valve_ring_pts(1,:) - ring_center(1); 
y = skeleton.valve_ring_pts(2,:) - ring_center(2); 

theta = atan2(y,x); 

[theta_unique idx] = unique(theta); 

ring_pts_unique = skeleton.valve_ring_pts(:,idx); 


% center mesh on zero, because atan2 returns results on this interval 
mesh = mod(mesh, 2*pi) - pi; 

theta_three_period = [theta_unique - 2*pi, theta_unique, theta_unique + 2*pi]; 

ring_points_three_period = [ring_pts_unique, ring_pts_unique, ring_pts_unique]; 

pts = interp1(theta_three_period, ring_points_three_period', mesh, 'spline')'; 
