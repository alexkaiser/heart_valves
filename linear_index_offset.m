function idx = linear_index_offset(j,k,N)
    %
    % Maps j,k in 3d to correct offset in flattened array 
    % Assumse N internal points with triangle domain 
    % 
      
    prev_values = N*(N+1)/2 - (N-k+1)*(N-k+2)/2; 
    idx = 3 * (j-1 + prev_values); 
end 


