function [is_internal is_bc linear_idx_offset point_idx_with_bc] = get_util_arrays_bead_slip(leaflet)
% 
% Returns three arrays with information about the geometry 
% 
% Output: 
%     is_internal          Boolean, true if 
%     is_bc                Point is a boundary condition that is fixed 
%     linear_idx_offset    In Jacobian, linear_idx_offset(j,k) + 1:3
%                          are the indices for the vector X(:,j,k)
% 

N                       = leaflet.N; 
j_max                   = leaflet.j_max; 
k_max                   = leaflet.k_max; 

is_internal       = zeros(j_max, k_max); 
is_bc             = zeros(j_max, k_max); 
linear_idx_offset = zeros(j_max, k_max); 
point_idx_with_bc = zeros(j_max, k_max); 


if leaflet.radial_and_circumferential 
    
    % radial and circumferential fiber layout 
    
    % valve ring at k_max
    k=k_max; 
    for j=1:j_max 
        is_bc(j,k) = true; 
    end 
    
    % loop from left free edge then up in k 
    j = k_max-1;  
    for k=1:(k_max-1)
        for k_tmp=k:(k_max-1)
            is_internal(j,k_tmp) = true; 
        end 
        j = j - 1; 
    end 

    
    j = k_max-1 + 1; 
    for k=1:(k_max-1)
        for k_tmp=k:(k_max-1)
            is_internal(j,k_tmp) = true; 
        end 
        j = j + 1; 
    end 
        
else 
    error('diag fibers not implemented with bead slip')
end 

count = 0; 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)
            linear_idx_offset(j,k) = count; 
            count = count + 3; 
        end 
    end 
end

count = 0;
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k) || is_bc(j,k)
            point_idx_with_bc(j,k) = count; 
            count = count + 1; 
        end 
    end 
end
