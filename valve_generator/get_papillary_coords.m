function papillary = get_papillary_coords(center, r, n_points, min_angle, max_angle)
% 
% Returns a distributed list of papillary points 
% 


angles = linspace(min_angle, max_angle, n_points); 
papillary = zeros(3,n_points); 

for i=1:n_points 
    papillary(:,i) = center + r * [cos(angles(i)); sin(angles(i)); 0]; 
end 


