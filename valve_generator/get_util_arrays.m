function [is_internal is_bc linear_idx_offset] = get_util_arrays(leaflet)
% 
% Returns three arrays with information about the geometry 
% 
% Output: 
%     is_internal          Boolean, true if 
%     is_bc                Point is a boundary condition that is fixed 
%     linear_idx_offset    In Jacobian, linear_idx_offset(j,k) + 1:3
%                          are the indices for the vector X(:,j,k)
% 

N = leaflet.N; 


if leaflet.radial_and_circumferential
    error('not implemented')
else 
    
    is_internal       = zeros(N+1, N+1); 
    is_bc             = zeros(N+1, N+1);
    linear_idx_offset = zeros(N+1, N+1);
    
    k = N+1; 
    for j=1:N+1
        is_bc(j,k) = true;
        k = k-1; 
    end 
    
    for j=1:N
        for k=1:N
            % in the triangle? 
            if ((j+k) < (N+2))
                is_internal(j,k) = true;  
            end 
        end 
    end 
    
    count = 0; 
    for k=1:N
        for j=1:N
            if is_internal(j,k)
                linear_idx_offset(j,k) = count; 
                count = count + 3; 
            end 
        end 
    end     
    
end 


