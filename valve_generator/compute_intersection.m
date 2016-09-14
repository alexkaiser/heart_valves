function coords = compute_intersection(X_flat, j, k, leaflet)
% 
% Computes the intersection of the rays which have coordinates j,k
% 
% Input: 
%     X_flat          2d surface including valve ring 
%     j,k             Indices to fill 
%     filter_params   paramters for the filter 
% 
% Output: 
%     coords          2d vector for intersection 
% 

a = leaflet.filter.a; 


% find the indices corresponding diagonal points which lie on the ray
j_left  = leaflet.N + 2 - k; 
k_right = leaflet.N + 2 - j; 

% intersection is a solution to a linear system 
A = [a + X_flat(1,j_left,k), a - X_flat(1,j,k_right); 
         X_flat(2,j_left,k),   - X_flat(2,j,k_right)]; 

rhs = [X_flat(1,j_left,k) - X_flat(1,j,k_right); 
       X_flat(2,j_left,k) - X_flat(2,j,k_right); ]; 

soln = A \ rhs; 
s = soln(1); 

coords = [-a;0]*s + (1 - s)*X_flat(:,j_left,k); 

t = soln(2); 
coords_alt = [a;0]*t + (1 - t)*X_flat(:,j,k_right);

tol = 1e2 * eps; 
if norm(coords - coords_alt) > tol 
    error('two formulas from line disagree'); 
end 



