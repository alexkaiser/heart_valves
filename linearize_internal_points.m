function v_linearized = linearize_internal_points(v, params)
%
%  Takes the internal values in X and arranges them in a linear array 
%  

% total internal points in triangular domain 
total_internal = 3*params.N*(params.N+1)/2;
idx = 1; 

v_linearized = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:params.N
    for j=1:params.N
        % in the triangle?
        if (j+k) < (params.N+2)
            v_linearized(idx + (0:2)) = v(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end 


