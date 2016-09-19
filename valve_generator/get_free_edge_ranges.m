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
    
    j_max = N; 
    k_max = N/2; 
    
    N = leaflet.N; 
    
    free_edge_idx_left  = zeros(N/2,   2); 
    free_edge_idx_right = zeros(N/2,   2); 
    chordae_idx_left    = zeros(N  , N/2); 
    chordae_idx_right   = zeros(N  , N/2); 
    
    j = N/2; 
    for k=1:(N/2)
        free_edge_idx_left(k,:) = [j; k]; 
        chordae_idx_left(j,k)   = k; 
        j = j - 1; 
    end 

    j = N/2 + 1; 
    for k=1:(N/2)
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


