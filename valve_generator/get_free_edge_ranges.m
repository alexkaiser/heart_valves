function [j_max k_max free_edge_idx_left free_edge_idx_right chordae_idx_left chordae_idx_right] = get_free_edge_ranges(leaflet)
% 
% Returns two 2d arrays of indices 
% 
% Input: 
%    leaflet    Parameter struct 
% 
% Output: 
%    free_edge_idx_left    2d array of left chordae indices of the implicitly defined leaves. 
%                          The i-th leaf is not included in the tree, 
%                          but instead on the leaflet as 
%                          X(:,free_edge_idx_left(:,1),free_edge_idx_left(:,1))  
%                          
%    free_edge_idx_right   2d array of right chordae indices  
%
%    chordae_idx_left      If chordae_idx_left(j,k) is nonzero, then this contains 
%                          the leaf index which is connected to X(:,j,k)
% 
 

if leaflet.radial_and_circumferential

    N = leaflet.N; 
    trapezoidal_flat_points = leaflet.trapezoidal_flat_points; 
    
    j_max = N + trapezoidal_flat_points; 
    k_max = N/2; % always this  
        
    free_edge_idx_left  = zeros(k_max, 2); 
    free_edge_idx_right = zeros(k_max, 2); 
    chordae_idx_left    = zeros(j_max, k_max); 
    chordae_idx_right   = zeros(j_max, k_max); 
    
    % Left free edge starts at (k_max,1)
    % and ends at (1,k_max) 
    j = k_max;  
    for k=1:k_max
        free_edge_idx_left(k,:) = [j; k]; 
        chordae_idx_left(j,k)   = k; 
        j = j - 1; 
    end 

    % Right free edge starts 1 + trapezoidal_flat_points 
    % to the right of left free edge corner 
    % and ends at (j_max, k_max)
    j = k_max + 1 + trapezoidal_flat_points; 
    for k=1:k_max
        free_edge_idx_right(k,:) = [j; k]; 
        chordae_idx_right(j,k)   = k; 
        j = j + 1; 
    end 

else 
    
    N = leaflet.N; 
    free_edge_idx_left  = [ones(N,1), (1:N)'];
    free_edge_idx_right = [(1:N)', ones(N,1)];

    j_max = N+1; 
    k_max = N+1; 
    
    chordae_idx_left  = zeros(j_max,k_max); 

    j=1; 
    for k=1:N
        chordae_idx_left(j,k) = k; 
    end 

    chordae_idx_right = zeros(N+1,N+1);

    k=1; 
    for j=1:N
        chordae_idx_right(j,k) = j; 
    end 

end 


