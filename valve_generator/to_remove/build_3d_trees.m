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










