function params = internal_points_to_2d(v_linearized, params)
%
%  Takes the internal values in X which are arranged in linear order
%  And places them back in the 3d vector array in params
%  

idx = 1; 

% here k is required to be the outer loop 
for k=1:params.N
    for j=1:params.N
        % in the triangle?
        if (j+k) < (N+2)
            params.X(:,j,k) = v_linearized(idx + (0:3)); 
            idx = idx + 3; 
        end 
    end 
end 




