function X_linearized = linearize_internal_points(leaflet, X, chordae)
%
%  Takes the internal values in X and arranges them in a linear array 
%  
%  If have a nonempty chordae data structure, then two additional arrays must be included
% 
%  Input: 
%      leaflet           Data parameters
%      v                 Three dimensional array 
%                        Has dimensions of leaflet
%                        Includes b.c.s and out of range data 
%      v_left_chordae    Left chordae tree if desired
%      v_right_chordae   Right chordae tree if desired
% 
%  Output: 
%      v_linearized      Internal points in a one dimensional array
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


% total internal points in triangular domain 

j_max       = leaflet.j_max; 
k_max       = leaflet.k_max; 
is_internal = leaflet.is_internal; 

total_internal = 3*sum(is_internal(:));
idx = 1; 

if (~exist('X', 'var')) && (~exist('X_chordae', 'var'))
    X       = leaflet.X; 
    chordae = leaflet.chordae; 
end 

if nargin == 2
    error('Must pass X and chordae or use current leaflet values'); 
end 

X_linearized = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:k_max
    for j=1:j_max
        if leaflet.is_internal(j,k)
            X_linearized(idx + (0:2)) = X(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end

% fprintf('length(X_learized) = %d, internal placed only\n', length(X_linearized)); 

for tree_idx = 1:leaflet.num_trees
    
    if isfield(leaflet, 'targets_for_bcs') && leaflet.targets_for_bcs         
        X_linearized = [X_linearized; chordae(tree_idx).root(:)];
        
        % fprintf('length(X_learized) = %d, tree_idx = %d, root placed\n', length(X_linearized), tree_idx); 
        
    end 
    
    X_linearized = [X_linearized; chordae(tree_idx).C(:)]; 
    % fprintf('length(X_learized) = %d, tree_idx = %d, tree placed\n', length(X_linearized), tree_idx); 
end 

% 'pause'

