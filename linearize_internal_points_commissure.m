function v_linearized = linearize_internal_points_commissure(v, params)
%
%  Takes the internal values in X and arranges them in a linear array 
%  

% total internal points in triangular domain 
total_internal = total_internal_commissure(params.N); 
idx = 1; 

v_linearized = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:params.N
    for j=1:params.N
        % in the triangle?
        if is_internal_commissure(j,k,params.N) 
            v_linearized(idx + (0:2)) = v(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end 


