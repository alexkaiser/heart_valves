function [j_max k_max free_edge_idx_left free_edge_idx_right chordae_idx_left chordae_idx_right] = get_free_edge_ranges_bead_slip(leaflet)
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
    k_max = N/2 + 1; 
        
    free_edge_idx_left  = zeros(N/2, 2); 
    free_edge_idx_right = zeros(N/2, 2); 
    chordae_idx_left    = zeros(j_max, k_max); 
    chordae_idx_right   = zeros(j_max, k_max); 
    
    % Left free edge starts at (N/2,1)
    % and ends at (1,N/2) 
    j = k_max-1;  
    for k=1:(k_max-1)
        free_edge_idx_left(k,:) = [j; k]; 
        chordae_idx_left(j,k)   = k; 
        j = j - 1; 
    end 

    % Right free edge starts N/2 + 1 
    % to the right of left free edge corner 
    % and ends at (j_max, k_max)
    j = (k_max-1) + 1; 
    for k=1:(k_max-1)
        free_edge_idx_right(k,:) = [j; k]; 
        chordae_idx_right(j,k)   = k; 
        j = j + 1; 
    end 

else 
    error('diagional not implemented for bead slip'); 
end 


