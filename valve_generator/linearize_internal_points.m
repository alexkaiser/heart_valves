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
N = leaflet.N; 
total_internal = 3*N*(N+1)/2;
idx = 1; 

v_linearized = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:N
    for j=1:N
        if leaflet.is_internal(j,k)
            v_linearized(idx + (0:2)) = v(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end

if exist('v_left_chordae', 'var') && exist('v_right_chordae', 'var')
    v_linearized = [v_linearized; v_left_chordae(:); v_right_chordae(:)]; 
end 