function v_linearized = linearize_internal_points_bead_slip_attached(valve, v_anterior, v_posterior, v_chordae_left, v_chordae_right)
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

j_max       = valve.anterior.j_max; 
k_max       = valve.anterior.k_max; 
is_internal = valve.anterior.is_internal; 

total_internal = 3*sum(is_internal(:));
idx = 1; 

v_linearized_anterior = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)
            v_linearized_anterior(idx + (0:2)) = v_anterior(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end


j_max       = valve.posterior.j_max; 
k_max       = valve.posterior.k_max; 
is_internal = valve.posterior.is_internal; 

total_internal = 3*sum(is_internal(:));
idx = 1; 

v_linearized_posterior = zeros(total_internal,1); 

% here k is required to be the outer loop 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)
            v_linearized_posterior(idx + (0:2)) = v_posterior(:,j,k); 
            idx = idx + 3; 
        end 
    end 
end

v_linearized = [v_linearized_anterior; v_linearized_posterior; v_chordae_left(:); v_chordae_right(:)]; 
