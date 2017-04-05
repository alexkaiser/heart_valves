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

for tree_idx = 1:leaflet.num_trees
    X_linearized = [X_linearized; chordae(tree_idx).C(:)]; 
end 

