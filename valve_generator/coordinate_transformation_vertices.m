function [vertices normal centroid R] = coordinate_transformation_vertices(vertices, ring_boundary_file, R_0, T_0)
% 
% takes file containing ring points 
% applies rigid translations and rotations that takes the unit circle to 
% the centroid and normal to the centroid and normal of the ring 

tol = 1e-14; 

[~, n_vertices] = size(vertices); 


f = fopen(ring_boundary_file, 'r'); 
vertices_ring_bdry = fscanf(f, '%f'); 
fclose(f); 

% first is number of vertices 
n_pts_ring_from_file = vertices_ring_bdry(1); 

% crop to keep just this 
vertices_ring_bdry = vertices_ring_bdry(2:end); 

vertices_ring_bdry = reshape(vertices_ring_bdry, 3, []); 
[~, n_pts_ring] = size(vertices_ring_bdry); 

if n_pts_ring ~= n_pts_ring_from_file
    error('inconsistent number of points from header and array size'); 
end 

centroid = zeros(3,1); 
for i=1:3
    centroid(i) = mean(vertices_ring_bdry(i,:)); 
end 

normal = cross(vertices_ring_bdry(:,1) - centroid, vertices_ring_bdry(:, floor(n_pts_ring/4)) - centroid); 
% take the z direction up normal 
if normal(3) < 0
    normal = -normal; 
end 
normal = normal / norm(normal); 


% initial rotation if requested 
if exist('T_0', 'var')
    vertices = T_0 + vertices; 
end 

% initial rotation if requested 
if exist('R_0', 'var')
    vertices = R_0 * vertices; 
end 

% check inverse rotation first 
% rotate normal so it has zero y component 
phi = atan2(normal(2), normal(3)); 
R_x = rotation_matrix_x(phi); 

normal_no_y = R_x * normal;  

if abs(normal_no_y(2)) > tol 
    error('did not remove y component'); 
end 

% then again so it has zero x 
theta = atan2(normal_no_y(1), normal_no_y(3)); 
R_y = rotation_matrix_y(theta); 

normal_zdir_from_rotation = R_y * normal_no_y; 

initial_normal = [0; 0; 1]; 

if norm(normal_zdir_from_rotation - initial_normal) > tol 
    error('inverse rotation incorrect'); 
end 

% rotation to apply is the inverse of these 
% and note that all matrices are orthogonal 
R = (R_x') * (R_y');   

if norm(R*initial_normal - normal) > tol
    error('R*[0;0;1] does not give correct normal'); 
end 


if norm(R*R' - eye(3)) > tol
    error('rotation matrix is not orthogonal'); 
end 

% apply rotation matrix to each vertex 
vertices = R * vertices; 

% place origin where desired 
for j=1:n_vertices
    vertices(:,j) = vertices(:,j) + centroid; 
end 
