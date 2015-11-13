function idx = linear_index_offset_commissure(j,k,N)
%
% Maps j,k in 3d to correct offset in flattened array 
% Assumse N internal points with triangle domain 
% 
      
if mod(N,2) ~= 1
    error('Must use odd N on commissural leaflets'); 
end 

% values in previous rows 
if k == 1
    prev_rows = 0; 
else 
    prev_rows = ((N+1)/2)^2 - ((N - 2*k + 3)/2)^2; 
end 

prev_this_row = j - k - 1; 

idx = 3*(prev_rows + prev_this_row); 



