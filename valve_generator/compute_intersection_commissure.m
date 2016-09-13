function coords = compute_intersection_commissure(X_flat, j, k, filter_params, left)
% 
% Computes the intersection of the rays which have coordinates j,k
% 
% Input: 
%     X_flat          2d surface including valve ring 
%     j,k             Indices to fill 
%     filter_params   paramters for the filter 
%     left            Uses left papillary if true
% 
% Output: 
%     coords          2d vector for intersection 
% 


if left 
    papillary = [-filter_params.a; 0]; 
else 
    papillary = [ filter_params.a; 0]; 
end 

N = filter_params.N; 
if mod(N,2) ~= 1
    error('must use odd N for commisural leaflet')
end

% k_vert is the k on the b.c. above the current j 
if j <= ((N+3)/2)
    k_vert = j; 
else 
    k_vert = N + 3 - j; 
end 

% j_horiz is the j on the right boundary 
j_horiz = N + 3 - k; 

% the left boundary is always k,k so no extra term needed 

% intersection is a solution to a linear system 
A = [X_flat(1,k,k) - X_flat(1,j_horiz,k), -papillary(1) + X_flat(1,j,k_vert);  
     X_flat(2,k,k) - X_flat(2,j_horiz,k), -papillary(2) + X_flat(2,j,k_vert)];   

rhs = [-X_flat(1,j_horiz,k) + X_flat(1,j,k_vert); 
       -X_flat(2,j_horiz,k) + X_flat(2,j,k_vert)]; 

soln = A \ rhs; 
s = soln(1); 

coords = s*X_flat(:,k,k) + (1 - s)*X_flat(:,j_horiz,k); 

t = soln(2); 
coords_alt = t*papillary + (1 - t)*X_flat(:,j,k_vert);

tol = 1e2 * eps; 
if norm(coords - coords_alt) > tol 
    error('two formulas from line disagree'); 
end 



