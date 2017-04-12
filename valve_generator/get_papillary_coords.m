function points = get_papillary_coords(valve, papillary_idx, n_points, min_angle, max_angle)
% 
% Returns a distributed list of papillary points 
% 


papillary = valve.skeleton.papillary(papillary_idx); 
c         = papillary.center; 
r         = papillary.radius; 
n         = papillary.normal; 

% this points from left to right papillary
% this is approximately in the 'y' direction, but determined by tip location 
approx_vertical = valve.skeleton.l_to_r_papillary; 

% project onto orthogonal complement of span(n)
% and normalize 
% local 'y' coord in the circle 
approx_vertical = (eye(3) - n*n') * approx_vertical; 
approx_vertical = approx_vertical / norm(approx_vertical); 

% local 'x' coord in the circle 
approx_horiz    = cross(approx_vertical, n); 

if abs(norm(approx_horiz) - 1) > eps 
    error('Should have gotten a unit vector by construction');     
end 

angles = linspace(min_angle, max_angle, n_points); 
points = zeros(3,n_points); 

for i=1:n_points 
    points(:,i) = c + ...
                  r * cos(angles(i)) * approx_horiz + ...
                  r * sin(angles(i)) * approx_vertical; 
end 





