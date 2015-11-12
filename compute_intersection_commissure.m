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

% intersection is a solution to a linear system 
A = [X_flat(1,k,k) - X_flat(1,N+3-k,k), -papillary(1) + X_flat(1,j,j);  
     X_flat(2,k,k) - X_flat(2,N+3-k,k), -papillary(2) + X_flat(2,j,j)];   

rhs = [-X_flat(1,N+3-k,k) + X_flat(1,j,j); 
       -X_flat(2,N+3-k,k) + X_flat(2,j,j)]; 

soln = A \ rhs; 
s = soln(1); 

coords = s*X_flat(:,k,k) + (1 - s)*X_flat(:,N+3-k,k); 

t = soln(2); 
coords_alt = t*papillary + (1 - t)*X_flat(:,j,j);

tol = 1e2 * eps; 
if norm(coords - coords_alt) > tol 
    error('two formulas from line disagree'); 
end 



