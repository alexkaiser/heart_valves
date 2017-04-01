function chordae = add_chordae(leaflet)
    % 
    % Adds a chordae data structure to the current parameters 
    % 
    % Input: 
    %     leaflet     Current leaflet data structures 
    % 
    % Output: 
    %     chordae     Chrodae tree data structure            
    % 

    % unpack relevant parameters 
    tree_frac           = leaflet.tree_frac; 
    left_papillary      = leaflet.left_papillary;
    right_papillary     = leaflet.right_papillary;
    X                   = leaflet.X; 
    k_multiplier        = leaflet.k_multiplier; 
    k_0                 = leaflet.k_0; 
    
    % free edge of leaflet at minimum k 
    k_min               = leaflet.k_min; 
    
    [n_leaves, m] = size(leaflet.free_edge_idx_left);  

    if m ~= 2
        error('free edge indices must be an N by 2 array'); 
    end 
        
    n_tree = log2(n_leaves);

    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 
    
%     if n_tree < 2
%         error('weird boundary errors possible on such a small tree'); 
%     end 
    
    if ~((0 < tree_frac) && (tree_frac < 1))
        error('multiplier on tree position must be between zero and one'); 
    end 

    
    total_len = 2^(n_tree+1) - 1; 
    max_internal = 2^(n_tree) - 1;     % last node that is not a leaf 
    
    C_left  = zeros(3,total_len); 
    C_right = zeros(3,total_len); 
        
    % initialize the left boundary conditions from the leaflet 
    for j=1:n_leaves 
        k_left = k_min(j); 
        C_left (:, max_internal + j)  = X(:,j,k_left);
        
        k_right = k_min(j + n_leaves); 
        C_right(:, max_internal + j)  = X(:, j + n_leaves, k_right); 
    end 
    
    for i=1:max_internal

        p = get_parent(C_left, i, left_papillary); 
        
        if ~is_leaf(i, max_internal)
            l = get_left__descendant(C_left, i, max_internal); 
            r = get_right_descendant(C_left, i, max_internal); 
            C_left(:,i) = (tree_frac)*p + 0.5*(1-tree_frac)*l + 0.5*(1-tree_frac)*r; 
        end
        
        p = get_parent(C_right, i, right_papillary); 
        
        if ~is_leaf(i, max_internal)
            l = get_left__descendant(C_right, i, max_internal); 
            r = get_right_descendant(C_right, i, max_internal); 
            C_right(:,i) = (tree_frac)*p + 0.5*(1-tree_frac)*l + 0.5*(1-tree_frac)*r;  
        end 
        
    end 
    
    % do not actually want the leaves (which are copied)
    % remove them 
    C_left  = C_left (:,1:max_internal); 
    C_right = C_right(:,1:max_internal); 
    
    
    % there are max_internal points on the tree 
    % leaves are not included as they are part of the leaflet 
    [m max_internal] = size(C_left); 
    
    % this parameter is the same N as in the leaflet 
    % it is the number of (not included) leaves in the tree
    NN = max_internal + 1; 
    
    % sanity check in building a balanced tree 
    n_tree = log2(NN);
    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 
        
    % each keeps parent wise spring constants 
    k_l = zeros(max_internal,1); 
    k_r = zeros(max_internal,1); 
        
    % set spring constant data structures 
    num_at_level = n_leaves/2; 
    idx = max_internal; 
    
    % constants connecting to the leaflet are inherited
    % first internal constant to the tree is k_multiplier times that 
    k_running = k_multiplier*k_0; 
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
    
    % reference configuration (ignoring the length multiplier) is the initial configuration 
    % copy it 
    chordae.Ref_l            = C_left; 
    chordae.Ref_r            = C_right; 
    
    chordae.k_l              = k_l; 
    chordae.k_r              = k_r; 
    chordae.k_0              = k_0; 
    chordae.k_multiplier     = k_multiplier; 
end 




function leaf = is_leaf(i, max_internal)
% checks whether node at indenx i is leaf 
% 
    if i > max_internal
        leaf = true; 
    else 
        leaf = false; 
    end 
end 


function r = get_right_descendant(C, i, max_internal)
% 
% finds the right most leaf which is a descendant of the node i 
% 
     
    if is_leaf(i, max_internal)
        error('trying to get descendants on a leaf'); 
    end 
    
    j = 2*i + 1;
    while ~is_leaf(j, max_internal)
        j = 2*j + 1; 
    end 
    
    r = C(:,j); 
end 


function l = get_left__descendant(C, i, max_internal)
% 
% finds the left most leaf which is a descendant of the node i 
% 
     
    if is_leaf(i, max_internal)
        error('trying to get descendants on a leaf'); 
    end 
    
    j = 2*i;
    while ~is_leaf(j, max_internal)
        j = 2*j; 
    end 
    
    l = C(:,j);
end 








