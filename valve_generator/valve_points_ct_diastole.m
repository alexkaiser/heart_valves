function skeleton = valve_points_ct_diastole(output)


% 
% Takes a center, radius and papillary coords 
% Moves center to origin 
% Takes papillary tips in rigid manner
% 

if ~exist('output', 'var')
    output = false; 
end 


radius = 0.1 * 16.065878777687722; 

normal = [-0.70870559417192203; 0.61353537460948204; 0.34829689187850299]; 


% initial coordinates are in mm
% multiply by .1 everywhere to get cm 

center = 0.1 * [43.211190372043184; 12.739459457869332; -18.011383112933899]; 
papillary_left  = 0.1 * [62.051506; 1.486666; -45.702660]; 
papillary_right = 0.1 * [84.484848; 3.326258; -21.130402]; 


% translate center to origin 
papillary_left  = papillary_left  - center; 
papillary_right = papillary_right - center; 

% get angle to rotate normal to in the y,z plane 
% and rotate accordingly 
theta = atan2(normal(2), normal(1)); 
R = rotation_matrix_z(-theta); 

normal          = R*normal; 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 

if abs(normal(2)) > eps
    error('did not land in the xz plane'); 
end 

% get angle to place on positive z axis 
phi = atan2(normal(3), normal(1)); 
R = rotation_matrix_y(pi/2 - phi); 

normal          = R*normal; 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 

if (abs(normal(1)) > eps) || (abs(normal(2)) > eps) 
    error('did not land in the z axis'); 
end 

if abs(normal(3) - 1) > eps
    error('z component of end normal is not in the right place'); 
end 

normal_to_papillary_plane = cross(papillary_left, papillary_right); 

% compute the midpoint
% rotate this to the negative x axis 
midpoint = 0.5 * (papillary_left(1:2) + papillary_right(1:2)); 
theta = atan2(midpoint(2), midpoint(1)); 

R = rotation_matrix_z(-theta + pi); 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 

radius 
papillary_left
papillary_right


th = 0:.1:2*pi; 
figure; 
plot3(radius*cos(th), radius*sin(th), zeros(size(th))); 
hold on 
plot3(papillary_left(1), papillary_left(2), papillary_left(3), '*'); 
plot3(papillary_right(1), papillary_right(2), papillary_right(3), 'o'); 
legend('ring', 'left', 'right'); 
title('mitral skeleton geometry'); 
xlabel('x')
ylabel('y')
zlabel('z')




skeleton.r                   = radius; 

tip_radius = .25; 

skeleton.papillary(1).name   = 'left'; 
skeleton.papillary(1).center = papillary_left;
skeleton.papillary(1).radius = left_radius;
skeleton.papillary(1).normal = left_normal;

skeleton.l_to_r_papillary    = l_to_r_papillary; 

skeleton.papillary(2).name   = 'right'; 
skeleton.papillary(2).center = papillary_right;
skeleton.papillary(2).radius = right_radius;
skeleton.papillary(2).normal = right_normal;





