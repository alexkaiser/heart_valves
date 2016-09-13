function v_linearized = linearize_internal_points(v, params, v_left_chordae, v_right_chordae)
%
%  Takes the internal values in X and arranges them in a linear array 
%  
%  If have a nonempty chordae data structure, then two additional arrays must be included
% 
%  Input: 
%      v                 Three dimensional array 
%                        Has dimensions of leaflet
%                        Includes b.c.s and out of range data 
%      params            Data parameters
%      v_left_chordae    Left chordae tree if desired
%      v_right_chordae   Right chordae tree if desired
% 
%  Output: 
%      v_linearized      Internal points in a one dimensional array
% 

% total internal points in triangular domain 
N = params.N; 
total_internal = 3*N*(N+1)/2;
idx = 1; 

v_linearized = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:N
    for j=1:N
        % in the triangle?
        if (j+k) < (N+2)
            v_linearized(idx + (0:2)) = v(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end

    
if isfield(params, 'chordae') 
    if ~isempty(params.chordae)
        
        if nargin < 4
            error('must include chordae arrays if they are to be linearized'); 
        end 
        
        v_linearized = [v_linearized; v_left_chordae(:); v_right_chordae(:)]; 
    end 
end 


