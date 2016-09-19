function coords = compute_point_of_leaflet(X_flat, leaflet)
% 
% Computes the intersection of the rays which have coordinates j,k
% Works for diagonal fibers only 
% 
% Input: 
%     X_flat          2d surface including valve ring 
%     j,k             Indices to fill 
%     filter_params   paramters for the filter 
% 
% Output: 
%     coords          2d vector for intersection 
% 

if ~leaflet.radial_and_circumferential
    error('Compute point assumes radial and circumferential geometry'); 
end 

a     = leaflet.filter.a; 
j_max = leaflet.j_max; 
k_max = leaflet.k_max; 


left_point  = X_flat(:,1,     k_max); 
right_point = X_flat(:,j_max, k_max); 

% intersection is a solution to a linear system 
A = [a + right_point(1), a - left_point(1); 
         right_point(2),   - left_point(2)]; 
     
rhs = [right_point(1) - left_point(1); 
       right_point(2) - left_point(2); ]; 


   
% % intersection is a solution to a linear system 
% A = [a + X_flat(1,j_left,k), a - X_flat(1,j,k_right); 
%          X_flat(2,j_left,k),   - X_flat(2,j,k_right)]; 
% 
% rhs = [X_flat(1,j_left,k) - X_flat(1,j,k_right); 
%        X_flat(2,j_left,k) - X_flat(2,j,k_right); ]; 

soln = A \ rhs; 
s = soln(1); 

coords = [-a;0]*s + (1 - s)*right_point; 

t = soln(2); 
coords_alt = [a;0]*t + (1 - t)*left_point;

tol = 1e2 * eps; 
if norm(coords - coords_alt) > tol 
    error('two formulas from line disagree'); 
end 

 