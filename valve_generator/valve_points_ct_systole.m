% 
% Takes a center, radius and papillary coords 
% Moves center to origin 
% Takes papillary tips in rigid manner
% 

radius = 0.1 * 21.885241; 


ring_point_one = 0.1 * [34.324305; -7.449943; -10.627514];
ring_point_two = 0.1 * [40.091951; 9.874523; -33.887690]; 


% initial coordinates are in mm
% multiply by .1 everywhere to get cm 

center          = 0.1 * [51.119208;  5.729583; -15.443644]; 
papillary_left  = 0.1 * [66.194740;  0.831810; -50.110394]; 
papillary_right = 0.1 * [85.581177; -1.730428; -27.283619]; 


normal = -cross(ring_point_one - center, ring_point_two - center);
normal = normal / norm(normal); 

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
axis equal 
hold on 
plot3(papillary_left(1), papillary_left(2), papillary_left(3), '*'); 
plot3(papillary_right(1), papillary_right(2), papillary_right(3), 'o'); 
legend('ring', 'left', 'right'); 
title('mitral skeleton geometry'); 
xlabel('x')
ylabel('y')
zlabel('z')






