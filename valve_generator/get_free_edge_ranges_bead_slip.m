function leaflet = get_free_edge_ranges_bead_slip(leaflet)
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
        
    free_edge_idx_left  = zeros(N/2, 2); 
    free_edge_idx_right = zeros(N/2, 2); 
    chordae_idx_left    = zeros(N, N/2 + 1); 
    chordae_idx_right   = zeros(N, N/2 + 1); 
    
    % Left free edge starts at (N/2,1)
    % and ends at (1,N/2) 
    j = k_max;  
    for k=1:k_max
        free_edge_idx_left(k,:) = [j; k]; 
        chordae_idx_left(j,k)   = k; 
        j = j - 1; 
    end 

    % Right free edge starts N/2 + 1 
    % to the right of left free edge corner 
    % and ends at (j_max, k_max)
    j = k_max + 1; 
    for k=1:k_max
        free_edge_idx_right(k,:) = [j; k]; 
        chordae_idx_right(j,k)   = k; 
        j = j + 1; 
    end 

else 
    error('diagional not implemented for bead slip'); 
end 



if isfield(leaflet, 'N_ring_to_ring') && leaflet.N_ring_to_ring > 0
    
    N_ring_to_ring = leaflet.N_ring_to_ring; 
    
    if N_ring_to_ring > (N/2 - 1)
        error('Not enough ring points to add that many ring to ring fibers'); 
    end 
    
    % ring_k_idx(j) tells us that X(j, ring_k_idx(j)) is a ring point 
    ring_k_idx = zeros(j_max,1); 
    
    k_idx = k_max + 1;
    % indices go up while fibers are being laid down 
    for j=1:N_ring_to_ring
        ring_k_idx(j) = k_idx;
        k_idx = k_idx + 1; 
    end 
    
    % indices stay same here, no more connecting fibers 
    if N_ring_to_ring < (N/2)
        for j=(N_ring_to_ring + 1):(j_max - N_ring_to_ring)
            ring_k_idx(j) = k_idx; 
        end 
    end 
    
    % indices go back down towards the other commissure 
    for j=(j_max - N_ring_to_ring + 1):j_max
        k_idx = k_idx - 1; 
        ring_k_idx(j) = k_idx; 
    end
    
    % arrays are now resized  
    k_max = k_max + N_ring_to_ring + 1; 
    
    % resize and zero pad chordae arrays
    chordae_idx_left (j_max,k_max) = 0; 
    chordae_idx_right(j_max,k_max) = 0; 
    
else 
    % extra row gets no additional fibers but placed for ring only 
    k_max = k_max + 1; 
    ring_k_idx = k_max * ones(j_max,1); 
end 


leaflet.j_max               = j_max; 
leaflet.k_max               = k_max; 
leaflet.free_edge_idx_left  = free_edge_idx_left; 
leaflet.free_edge_idx_right = free_edge_idx_right; 
leaflet.chordae_idx_left    = chordae_idx_left; 
leaflet.chordae_idx_right   = chordae_idx_right; 
leaflet.ring_k_idx          = ring_k_idx; 


