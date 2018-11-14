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

% periodic wrap 
% x = [x, x(1)]; 
% y = [y, y(1)]; 

theta = atan2(y,x); 

% center mesh on zero, because atan2 returns results on this interval 
mesh = mod(mesh, 2*pi) - pi; 

pts_x = interp1(theta, skeleton.valve_ring_pts(1,:), mesh)'; 
pts_y = interp1(theta, skeleton.valve_ring_pts(2,:), mesh)'; 
pts_z = interp1(theta, skeleton.valve_ring_pts(3,:), mesh)'; 

pts = [pts_x, pts_y, pts_z]; 

