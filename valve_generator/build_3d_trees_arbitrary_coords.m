function chordae = build_3d_trees_arbitrary_coords(params, filter_params, tree_frac, k_0, k_multiplier) 

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
    
    
  % 
    % Creates and returns two binary trees of chordae 
    % The binary trees are fully balanced by definition
    % N, the number of internal nodes, must be a power of two  
    % 
    % Input 
    % params                  
    % left_papillary_flat     Left papillary muscle location 
    % right_papillary_flat    Right papillary muscle location 
    % tree_frac               Coefficient of parent in weighted sum for placement of parents
    %                          (tree_frac)*p + 0.5*(1-tree_frac)*l + 0.5*(1-tree_frac)*r
    % 
    

    % Copyright (c) 2019, Alexander D. Kaiser
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are met:
    % 
    % 1. Redistributions of source code must retain the above copyright notice, this
    %    list of conditions and the following disclaimer.
    % 
    % 2. Redistributions in binary form must reproduce the above copyright notice,
    %    this list of conditions and the following disclaimer in the documentation
    %    and/or other materials provided with the distribution.
    % 
    % 3. Neither the name of the copyright holder nor the names of its
    %    contributors may be used to endorse or promote products derived from
    %    this software without specific prior written permission.
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


    [X,alpha,beta,N,p_0,R,ref_frac] = unpack_params(params);  
   
    n_tree = log2(N);
    n_leaves = N; 
    
    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 
    
    if n_tree < 2
        error('weird boundary errors possible on such a small tree'); 
    end 
    
    if ~((0 < tree_frac) && (tree_frac < 1))
        error('multiplier on tree position must be between zero and one'); 
    end 

    if isfield(filter_params, 'left_papillary') && (isfield(filter_params, 'right_papillary'))
        left_papillary  = filter_params.left_papillary;
        right_papillary = filter_params.right_papillary;
    else 
        error('Cannot call arbitrary coords tree builder without papillary coords set'); 
    end 
    
    
    total_len = 2^(n_tree+1) - 1; 
    max_internal = 2^(n_tree) - 1;     % last node that is not a leaf 
    
    C_left  = zeros(3,total_len); 
    C_right = zeros(3,total_len); 
        
    % initialize the left boundary conditions from the leaflet 
    for i=1:n_leaves
        C_left (:, max_internal + i)  = X(:,1,i); 
        C_right(:, max_internal + i)  = X(:,i,1); 
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
    

       
    % Reset size parameters, now in the 3d part...  
    % Fix/check this later
    
    
    
    
    % there are max_internal points on the tree 
    % leaves are not included as they are part of the leaflet 
    [m max_internal] = size(C_left); 
    
    % this parameter is the same N as in the leaflet 
    % it is the number of (not included) leaves in the tree
    N = max_internal + 1; 
    
    % sanity check in building a balanced tree 
    n_tree = log2(N);
    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 
        
    % each keeps parent wise spring constants 
    k_l = zeros(max_internal,1); 
    k_r = zeros(max_internal,1); 
        
    % set spring constant data structures 
    num_at_level = N/2; 
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






