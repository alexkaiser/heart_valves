function chordae = build_3d_trees(chordae_flat, filter_params, k_0, k_multiplier) 

    % 
    % Creates and returns two binary trees of chordae 
    % The binary trees are fully balanced by definition
    % N, the number of internal nodes, must be a power of two  
    % 
    % Input 
    % chordae_flat        Data structure of surface 
    % filter_params       Current filter parameters 
    % k_0                 Spring constant on leaves 
    % k_multiplier        Each additional level of the tree is multiplied by k_multiplier 
    % 
    

    
    [C_left_flat, C_right_flat, left_papillary_flat, right_papillary_flat] = unpack_chordae_flat(chordae_flat); 
    
    % there are max_internal points on the tree 
    % leaves are not included as they are part of the leaflet 
    [m max_internal] = size(C_left_flat); 
    
    if m ~= 2
        error('flat chordae must be two dimensional')
    end 
    
    % this parameter is the same N as in the leaflet 
    % it is the number of (not included) leaves in the tree
    N = max_internal + 1; 
    
    % sanity check in building a balanced tree 
    n_tree = log2(N);
    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 
        
    C_left = zeros(3,max_internal); 
    C_right = zeros(3,max_internal); 
    
    % each keeps parent wise spring constants 
    k_l = zeros(max_internal,1); 
    k_r = zeros(max_internal,1); 
        
    % evaluate the map on the 2d nodes 
    left_papillary  = cone_filter( left_papillary_flat(1),  left_papillary_flat(2), filter_params); 
    right_papillary = cone_filter(right_papillary_flat(1), right_papillary_flat(2), filter_params); 
    
    for i=1:max_internal
        C_left (:,i) = cone_filter( C_left_flat(1,i),  C_left_flat(2,i), filter_params); 
        C_right(:,i) = cone_filter(C_right_flat(1,i), C_right_flat(2,i), filter_params); 
    end 
        
    % set refernce lengths 
%     for i=1:total_len
%         p = get_parent(C_left, i, left_papillary);         
%         Ref_l(i) = len_multiplier * norm(C_left(:,i) - p);
%         
%         p = get_parent(C_right, i, right_papillary); 
%         Ref_r(i) = len_multiplier * norm(C_right(:,i) - p); 
%     end 
    
    % set spring constant data structures 
    num_at_level = N/2; 
    idx = max_internal; 
    
    % constants connecting to the leaflet are inherited
    % first internal constant to the tree is twice that 
    k_running = 2*k_0; 
    while num_at_level >= 1
    
        for j=1:num_at_level
            k_l(idx) = k_running; 
            k_r(idx) = k_running; 
            idx = idx - 1; 
        end 
        
        k_running = k_running * k_multiplier; 
        num_at_level = num_at_level / 2; 
    end 
    
    % set up return structure 
    chordae.C_left           = C_left; 
    chordae.C_right          = C_right; 
    chordae.left_papillary   = left_papillary;
    chordae.right_papillary  = right_papillary;
%     chordae.Ref_l            = Ref_l; 
%     chordae.Ref_r            = Ref_r; 
    chordae.k_l              = k_l; 
    chordae.k_r              = k_r; 
    chordae.k_0              = k_0; 
    chordae.k_multiplier     = k_multiplier; 
end 


% function p = get_parent(C, i, papillary)
% % returns parent coordinates
% % if i==1 then returns papillary 
% 
%     if i==1
%         p = papillary; 
%     else 
%         p = C(:,floor(i/2));  
%     end 
% end 
% 
% 









