function skeleton = valve_points_ct_systole(output)


% 
% Takes a center, radius and papillary coords 
% Moves center to origin 
% Takes papillary tips in rigid manner
% 


if ~exist('output', 'var')
    output = false; 
end 


radius = 0.1 * 21.885241; 


ring_point_one = 0.1 * [34.324305; -7.449943; -10.627514];
ring_point_two = 0.1 * [40.091951; 9.874523; -33.887690]; 


% initial coordinates are in mm
% multiply by .1 everywhere to get cm 

center          = 0.1 * [51.119208;  5.729583; -15.443644]; 
papillary_left  = 0.1 * [66.194740;  0.831810; -50.110394]; 
papillary_right = 0.1 * [85.581177; -1.730428; -27.283619]; 


% polygons from papillary segmentation 
left_points = 0.1 * [68.126724 0.328857 -47.740620
                     64.706276 -1.425585 -50.922058
                     65.751228 3.592152 -51.668495]'; 
           
right_points = 0.1 * [88.904839 -1.376166 -23.674393
                      86.042999 -5.822769 -25.194181
                      83.157478 -4.150777 -29.086617
                      83.295311 0.410899 -30.680325
                      86.505264 2.286659 -27.782587]'; 

                  
[m n] = size(left_points); 
if m ~= 3
    error('Must always have three d vectors'); 
end 

radii = zeros(n,1); 

for j=1:n
    radii(j) = norm(papillary_left - left_points(:,j)); 
end 

left_radius = mean(radii); 
                  
left_normal = cross(left_points(:,2) - papillary_left, left_points(:,1) - papillary_left); 
left_normal = left_normal / norm(left_normal); 

[m n] = size(right_points); 
if m ~= 3
    error('Must always have three d vectors'); 
end 

radii = zeros(n,1); 

for j=1:n
    radii(j) = norm(papillary_right - right_points(:,j)); 
end 

right_radius = mean(radii); 

right_normal = cross(right_points(:,2) - papillary_right, right_points(:,1) - papillary_right);                  

right_normal = right_normal / norm(right_normal); 
                  

% sign chosen manually, since do not know orientation of original wrap 
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
left_normal     = R*left_normal; 
right_normal    = R*right_normal; 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 

if abs(normal(2)) > eps
    error('did not land in the xz plane'); 
end 

% get angle to place on positive z axis 
phi = atan2(normal(3), normal(1)); 
R = rotation_matrix_y(pi/2 - phi); 

normal          = R*normal; 
left_normal     = R*left_normal; 
right_normal    = R*right_normal; 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 

if (abs(normal(1)) > eps) || (abs(normal(2)) > eps) 
    error('did not land in the z axis'); 
end 

if abs(normal(3) - 1) > eps
    error('z component of end normal is not in the right place'); 
end 


% compute the midpoint
% rotate this to the negative x axis 
midpoint = 0.5 * (papillary_left(1:2) + papillary_right(1:2)); 
theta = atan2(midpoint(2), midpoint(1)); 

R = rotation_matrix_z(-theta + pi); 
papillary_left  = R*papillary_left; 
papillary_right = R*papillary_right; 

if output 

    radius 
    papillary_left
    left_radius
    left_normal 
    papillary_right
    right_radius
    right_normal

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
    
end 


l_to_r_papillary = papillary_right - papillary_left; 
l_to_r_papillary = l_to_r_papillary / norm(l_to_r_papillary);



skeleton.r                   = radius; 

skeleton.papillary(1).name   = 'left'; 
skeleton.papillary(1).center = papillary_left;
skeleton.papillary(1).radius = left_radius;
skeleton.papillary(1).normal = left_normal;

skeleton.l_to_r_papillary    = l_to_r_papillary; 

skeleton.papillary(2).name   = 'right'; 
skeleton.papillary(2).center = papillary_right;
skeleton.papillary(2).radius = right_radius;
skeleton.papillary(2).normal = right_normal;


