function idx = vertex_index_offset(j,k,N)
%
% Maps j,k in 3d to correct offset in flattened array 
% Assumse N internal points with triangle domain 
% Includes boundary points as valid 
% 
% This subtracts one because IBAMR uses zero indexing 

prev_values = (N+1)*(N+2)/2 - (N-k+2)*(N-k+3)/2; 
idx = j + prev_values - 1; 

