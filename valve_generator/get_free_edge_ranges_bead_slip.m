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
 
if leaflet.num_trees ~= 2
    error('not implemented')
end 


if leaflet.radial_and_circumferential

    N = leaflet.N; 
    
    j_max = N; 
    
    % minimum index is variable 
    k_min = zeros(j_max,1); 
    
    k_max = N/2; 
        
    for tree_idx = 1:leaflet.num_trees
        chordae(1).free_edge_idx  = zeros(N/2, 2); 
        chordae(2).free_edge_idx  = zeros(N/2, 2); 
    end 
    
    chordae_idx_left    = zeros(N, N/2 + 1); 
    chordae_idx_right   = zeros(N, N/2 + 1); 
    
    % Left free edge starts at (N/2,1)
    % and ends at (1,N/2) 
    j = k_max;  
    for k=1:k_max
        chordae(1).free_edge_idx(k,:) = [j; k];
        chordae_idx_left(j,k)   = k; 
        k_min(j) = k; 
        j = j - 1; 
    end 

    % Right free edge starts N/2 + 1 
    % to the right of left free edge corner 
    % and ends at (j_max, k_max)
    j = k_max + 1; 
    for k=1:k_max
        chordae(2).free_edge_idx(k,:) = [j; k]; 
        chordae_idx_right(j,k)   = k; 
        k_min(j) = k; 
        j = j + 1; 
    end 

else 
    error('diagional not implemented for bead slip'); 
end 



if isfield(leaflet, 'ring_to_ring_range') && (~isempty(leaflet.ring_to_ring_range)) && (max(leaflet.ring_to_ring_range) > 0) 
    
    if length(leaflet.ring_to_ring_range) == 1  
        min_ring_to_ring = 1; 
        max_ring_to_ring = leaflet.ring_to_ring_range; 
    elseif length(leaflet.ring_to_ring_range) == 2  
        min_ring_to_ring = leaflet.ring_to_ring_range(1); 
        max_ring_to_ring = leaflet.ring_to_ring_range(2); 
    else 
        error('Must specify exactly one or two indices for ring to ring fibers'); 
    end 
    
    % number of fibers placed 
    % N_ring_to_ring = 
    
    if max_ring_to_ring > (N/2 - 1)
        error('Not enough ring points to add that many ring to ring fibers'); 
    end 
    
    % ring_k_idx(j) tells us that X(j, ring_k_idx(j)) is a ring point 
    ring_k_idx = zeros(j_max,1); 
    
    % first k_idx is one up from internal 
    k_idx = k_max + 1;
    
    % points on ring with no connecting fibers 
    for j=1:(min_ring_to_ring-1)
        ring_k_idx(j) = k_idx; 
    end 
    
    % indices go up while fibers are being laid down 
    for j=min_ring_to_ring:max_ring_to_ring
        ring_k_idx(j) = k_idx;
        k_idx = k_idx + 1; 
    end 
    
    % indices stay same here, no more connecting fibers 
    if max_ring_to_ring < (N/2)
        for j=(max_ring_to_ring + 1):(j_max - max_ring_to_ring)
            ring_k_idx(j) = k_idx; 
        end 
    end 
    
    % indices go back down towards the other commissure 
    for j=(j_max - max_ring_to_ring + 1):(j_max - min_ring_to_ring + 1)
        k_idx = k_idx - 1; 
        ring_k_idx(j) = k_idx; 
    end
    
    % points on ring with no connecting fibers at other end 
    for j=(j_max - min_ring_to_ring + 2):j_max
        ring_k_idx(j) = k_idx; 
    end 
    
    
    % arrays are now resized  
    k_max = k_max + (max_ring_to_ring - min_ring_to_ring) + 2; 
    
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
leaflet.k_min               = k_min; 
leaflet.chordae             = chordae; 
leaflet.chordae_idx_left    = chordae_idx_left; 
leaflet.chordae_idx_right   = chordae_idx_right; 
leaflet.ring_k_idx          = ring_k_idx; 


