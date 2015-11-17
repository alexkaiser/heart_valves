function v_linearized = linearize_internal_points_commissure(v, params)
%
%  Takes the internal values in X and arranges them in a linear array 
%  

% total internal points in triangular domain 
N = params.N; 
total_internal = total_internal_commissure(N); 
idx = 1; 

v_linearized = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:((N+3)/2)
    for j=1:N+2
        if is_internal_commissure(j,k,N)
            v_linearized(idx + (0:2)) = v(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end 


