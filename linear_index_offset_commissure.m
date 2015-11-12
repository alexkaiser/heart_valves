function idx = linear_index_offset_commissure(j,k,N)
%
% Maps j,k in 3d to correct offset in flattened array 
% Assumse N internal points with triangle domain 
% 
      
if mod(N,2) ~= 1
    error('Must use odd N on commissural leaflets'); 
end 

prev_values = ((N+1)/2)^2 - ((N - 2*k + 1)/2)^2 - k; 
idx = 3*prev_values; 



