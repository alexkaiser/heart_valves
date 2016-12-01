function v_linearized = linearize_internal_points(leaflet, v, v_left_chordae, v_right_chordae)
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

v_linearized = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:k_max
    for j=1:j_max
        if leaflet.is_internal(j,k)
            v_linearized(idx + (0:2)) = v(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end

if exist('v_left_chordae', 'var') && exist('v_right_chordae', 'var') && ~leaflet.leaflet_only
    v_linearized = [v_linearized; v_left_chordae(:); v_right_chordae(:)]; 
end 