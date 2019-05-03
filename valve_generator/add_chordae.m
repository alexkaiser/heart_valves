function leaflet = add_chordae(leaflet, tree_idx)
    % 
    % Adds a chordae data structure to the current parameters 
    % 
    % Input: 
    %     leaflet     Current leaflet data structures 
    % 
    % Output: 
    %     chordae     Chrodae tree data structure            
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


    % unpack relevant parameters
    papillary           = leaflet.papillary;
    X                   = leaflet.X; 
    chordae             = leaflet.chordae; 
    free_edge_idx       = chordae(tree_idx).free_edge_idx; 
    
    [n_leaves, m] = size(free_edge_idx);  
    
    if m ~= 2
        error('free edge indices must be an N by 2 array'); 
    end 
        
    n_tree = log2(n_leaves);

    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 
    
    if n_tree < 1
        warning('weird boundary errors possible on such a small tree'); 
    end 

    if isfield(leaflet, 'targets_for_bcs') && leaflet.targets_for_bcs 
        targets_for_bcs = true; 
    else 
        targets_for_bcs = false; 
    end 
    
    % just for initial guess 
    tree_frac           = 0.5; 
    
    total_len = 2^(n_tree+1) - 1; 
    max_internal = 2^(n_tree) - 1;     % last node that is not a leaf 
   
    if ~targets_for_bcs
        % root pulled from the papillary arrays 
        chordae(tree_idx).root = papillary(:,tree_idx); 
    else 
        % root target pulled from the papillary arrays 
        chordae(tree_idx).root_target = papillary(:,tree_idx); 
        chordae(tree_idx).targets_for_bcs = true; 
        
        % just move the initial condition up slightly for the root 
        % which is now an internal point 
        chordae(tree_idx).root        = papillary(:,tree_idx) + [0; 0; 0.1]; 
    end 
    
    % initialize the left boundary conditions from the leaflet     
    chordae(tree_idx).C = zeros(3,total_len);   

    % Set the free edge according to array of indicies 
    for i=1:n_leaves 
        j = free_edge_idx(i,1); 
        k = free_edge_idx(i,2); 
        chordae(tree_idx).C(:, max_internal + i)  = X(:,j,k);
    end 
     
    
    for i=1:max_internal

        p = get_parent(chordae(tree_idx).C, i, chordae(tree_idx).root); 

        if ~is_leaf(i, max_internal)
            l = get_left__descendant(chordae(tree_idx).C, i, max_internal); 
            r = get_right_descendant(chordae(tree_idx).C, i, max_internal); 
            chordae(tree_idx).C(:,i) = (tree_frac)*p + 0.5*(1-tree_frac)*l + 0.5*(1-tree_frac)*r; 
        end

    end 
    
    
    % do not actually want the leaves (which are copied)
    % remove them 
    chordae(tree_idx).C = chordae(tree_idx).C(:,1:max_internal); 
   
    % there are max_internal points on the tree 
    % leaves are not included as they are part of the leaflet 
    [m max_internal] = size(chordae(tree_idx).C); 
    
    % this parameter is the number of (not included) leaves in the tree
    NN = max_internal + 1; 
    
    % sanity check in building a balanced tree 
    n_tree = log2(NN);
    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 
    
    chordae(tree_idx).min_global_idx   = leaflet.total_internal_with_trees + 1; 
    
    leaflet.total_internal_with_trees  = leaflet.total_internal_with_trees + 3*max_internal; 
    
    % one extra internal here for root
    if targets_for_bcs
        leaflet.total_internal_with_trees = leaflet.total_internal_with_trees + 3; 
    end 
    
    leaflet.chordae = chordae;
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








