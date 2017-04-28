function [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr, k_nbr)
%
% Returns whether neighbor is a valid point, 
% If valid retuns neighbor indicies 
% and incides for spring constants 
%
% 

j_max       = leaflet.j_max; 
k_max       = leaflet.k_max; 
is_internal = leaflet.is_internal; 
is_bc       = leaflet.is_bc;

if isfield(leaflet, 'periodic_j')
    periodic_j = leaflet.periodic_j; 
else
    periodic_j = zeros(k_max,1); 
end

% j spring is minimum, unless zero in which case gets a periodic wrap 
j_spr = min(j, j_nbr); 
if j_spr == 0 
    j_spr = j_max; 
end 

% j_nbr may need periodic reduction 
j_nbr = get_j_nbr(j_nbr, k, periodic_j, j_max); 

% k_nbr is always an identity operation, no periodicity in this direction 
% only include for consistency 
% k_nbr = k_nbr; 

k_spr = min(k, k_nbr);

% neighbor must be valid 
if (j_nbr > 0) && (k_nbr > 0) && ...
   (j_nbr <= j_max) && (k_nbr <= k_max) && ...
   (is_internal(j_nbr,k_nbr) || is_bc(j_nbr,k_nbr))
    valid = true; 
else 
    valid = false; 
end 
