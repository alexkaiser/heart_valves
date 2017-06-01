function leaflet = get_free_edge_ranges_bead_slip(leaflet, tree_n_start)
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
N                 = leaflet.N; 


j_max = N; 

% minimum index is variable 
k_chordae_attachment = zeros(j_max,1); 

% max possible k
k_max = N/2; 

k = k_max + leaflet_N_start; 

j = 1; 

% follow up/down direction of to variable minimum boundary of leaflet mesh 
for i = 1:length(N_per_direction)

    N_tmp     = N_per_direction(i); 
    direction = leaflet_direction(i); 

    for tmp = 1:N_tmp

        k_chordae_attachment(j) = k; 

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

j = tree_n_start; 

for tree_idx = 1:num_trees 

    n_leaves_tmp  = n_leaves(tree_idx); 

    chordae(tree_idx).free_edge_idx = zeros(n_leaves_tmp,2); 

    for leaf_idx=1:n_leaves_tmp

        k = k_chordae_attachment(j); 

        chordae(tree_idx).free_edge_idx(leaf_idx,:) = [j; k];
        chordae_idx(j,k).tree_idx = tree_idx;  
        chordae_idx(j,k).leaf_idx = leaf_idx; 

        % Horizonal index always increases 
        j = j + 1; 

        % periodic reduction if necessary 
        if j == (j_max + 1)
            j = 1;
        end 

    end 

end 


k_max = k_max + 1 + n_rings_periodic; 
ring_k_idx = k_max * ones(j_max,1); 

% resize and zero pad chordae arrays
chordae_idx(j_max,k_max).tree_idx = 0; 
chordae_idx(j_max,k_max).leaf_idx = 0;

periodic_j = zeros(k_max, 1); 

% minimum point up to bc at idx 1 gets periodic connection 
% can set or ignore the commissure point that is part of the leaflet here 
leaflet.comm_point_attached = false; 

if leaflet.comm_point_attached 
    error('This messes with lots of other assumptions. Not implemented for now.'); 
    k_start_periodic = k_chordae_attachment(1); 
else 
    k_start_periodic = k_chordae_attachment(1) + 1;
end 

% here we add points below the chordae attachment free edge if requested
if leaflet.n_edge_connectors > 0
    k_start_periodic = k_start_periodic - leaflet.n_edge_connectors; 
end 

% minimum valid is minimum chordae attachment
% or minimum with edge connectors 
k_min = min(k_start_periodic, k_chordae_attachment); 

for k=k_start_periodic:k_max
    periodic_j(k) = 1; 
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


leaflet.j_max                = j_max; 
leaflet.k_max                = k_max; 
leaflet.k_min                = k_min; 
leaflet.k_chordae_attachment = k_chordae_attachment; 
leaflet.periodic_j           = periodic_j; 
leaflet.chordae              = chordae; 
leaflet.chordae_idx          = chordae_idx; 
leaflet.ring_k_idx           = ring_k_idx; 


