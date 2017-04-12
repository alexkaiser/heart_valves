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
 
num_trees         = leaflet.num_trees; 
n_leaves          = leaflet.n_leaves; 
leaflet_direction = leaflet.leaflet_direction; 
N_per_direction   = leaflet.N_per_direction; 
leaflet_N_start   = leaflet.leaflet_N_start; 

n_rings_periodic  = leaflet.n_rings_periodic; 


if leaflet.radial_and_circumferential

    N = leaflet.N; 
    
    j_max = N; 
    
    % minimum index is variable 
    k_min = zeros(j_max,1); 
    
    % max possible k
    k_max = N/2; 
    
    k = k_max + leaflet_N_start; 
    
    j = 1; 
    
    % follow up/down direction of to variable minimum boundary of leaflet mesh 
    for i = 1:length(N_per_direction)
        
        N_tmp     = N_per_direction(i); 
        direction = leaflet_direction(i); 
        
        for tmp = 1:N_tmp
        
            k_min(j) = k; 
        
            % j always incremented 
            j = j+1; 
            
            % Do not increment last k 
            if tmp < N_tmp
                k = k + direction; 
            end        
        
        end 
        
    end 
    
        
    for tree_idx = 1:num_trees
        chordae(tree_idx).free_edge_idx  = zeros(N/2, 2); 
    end 
    
    chordae_idx(N,N/2 + 1).tree_idx = 0;  
    chordae_idx(N,N/2 + 1).leaf_idx = 0;  
    
    j = 1; 
    
    for tree_idx = 1:num_trees 
        
        n_leaves_tmp  = n_leaves(tree_idx); 
        
        chordae(tree_idx).free_edge_idx = zeros(n_leaves_tmp,2); 
        
        for leaf_idx=1:n_leaves_tmp
            
            k = k_min(j); 
        
            chordae(tree_idx).free_edge_idx(leaf_idx,:) = [j; k];
            chordae_idx(j,k).tree_idx = tree_idx;  
            chordae_idx(j,k).leaf_idx = leaf_idx;

%             
%             k_min(j) = k; 
%             
%             % Incremented in direction of sign
%             % Except on final iteration 
%             if leaf_idx < n_leaves_tmp
%                 k = k + leaflet_direction(tree_idx); 
%             end 
            
            % Horizonal index always increases 
            j = j + 1; 
            
        end 
        
    end 

else 
    error('diagional not implemented for bead slip'); 
end 



if isfield(leaflet, 'ring_to_ring_range') && (~isempty(leaflet.ring_to_ring_range)) && (max(leaflet.ring_to_ring_range) > 0) 
    
    if n_rings_periodic ~= 0
        error('ring to ring and periodic loops not implemented at same time'); 
    end 
    
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
    chordae_idx(j_max,k_max).tree_idx = 0; 
    chordae_idx(j_max,k_max).leaf_idx = 0; 
    
else 
    % extra row gets no additional fibers but placed for ring only 
    
    k_max_internal = k_max; 
    
    k_max = k_max + 1 + n_rings_periodic; 
    ring_k_idx = k_max * ones(j_max,1); 
    
    % resize and zero pad chordae arrays
    chordae_idx(j_max,k_max).tree_idx = 0; 
    chordae_idx(j_max,k_max).leaf_idx = 0;
    
    periodic_j = zeros(k_max, 1); 
    
    % minimum point up to bc at idx 1 gets periodic connection 
    for k=k_min(1):k_max
        periodic_j(k) = 1; 
    end 
    
end 


% clean up empty values in tension struct 
for j=1:j_max
    for k=1:k_max
        if isempty(chordae_idx(j,k).tree_idx)
            chordae_idx(j,k).tree_idx = 0; 
        end
        
        if isempty(chordae_idx(j,k).leaf_idx)
            chordae_idx(j,k).leaf_idx = 0; 
        end 
    end 
end 


leaflet.j_max               = j_max; 
leaflet.k_max               = k_max; 
leaflet.k_min               = k_min; 
leaflet.periodic_j          = periodic_j; 
leaflet.chordae             = chordae; 
leaflet.chordae_idx         = chordae_idx; 
leaflet.ring_k_idx          = ring_k_idx; 


