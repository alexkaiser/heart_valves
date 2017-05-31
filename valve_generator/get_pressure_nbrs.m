function [j_plus__1 j_minus_1 k_plus__1 k_minus_1 m] = get_pressure_nbrs(leaflet,j,k)
% 
% Returns neighbors for computing pressure  
% The point j,k should be internal (not a boundary condition)
% 
% If the point is not on the free edge, then points up and down
% for a two mesh widths wide, centered finite difference are supplied 
% 
% If the point is on the free edge, then a 
% 
% The value "m" is the coefficient of the finite difference pressure 
% This is four internal to the leaflet
% but is 1/2 if one direction uses a one mesh width stencil  
% and 1 if both directions use one mesh width stencils 
% 

num_j_modified = 0; 
[valid j_plus__1] = get_indices(leaflet, j, k, j+1, k); 
if ~valid
    j_plus__1 = j; 
    num_j_modified = num_j_modified + 1; 
end 

[valid j_minus_1] = get_indices(leaflet, j, k, j-1, k); 

if ~valid
    j_minus_1 = j; 
    num_j_modified = num_j_modified + 1; 
end 

num_k_modified = 0; 
[valid tmp k_plus__1] = get_indices(leaflet, j, k, j, k+1); 
if ~valid
    k_plus__1 = k; 
    num_k_modified = num_k_modified + 1; 
end 

[valid tmp k_minus_1] = get_indices(leaflet, j, k, j, k-1); 

if ~valid
    k_minus_1 = k; 
    num_k_modified = num_k_modified + 1; 
end 


if     (num_j_modified == 0) && (num_k_modified == 0)
    m = 0.25; 
elseif (num_j_modified == 1) && (num_k_modified == 0)
    m = 0.5; 
elseif (num_j_modified == 0) && (num_k_modified == 1)
    m = 0.5; 
elseif (num_j_modified == 1) && (num_k_modified == 1)
    m = 1.0; 
else 
    error('All pressure indices must be modified zero or one'); 
end 




