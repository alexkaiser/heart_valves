function J = build_jacobian(params, filter_params)
% 
% Builds the Jacobian for the current index and parameter values 
% 
% Input 
%      params    Current parameter values
%      j,k       Indices, must be internal to the arrays
% 
% Output 
%      J         Jacobian of difference equations 


[X,alpha,beta,N,p_0,R,ref_frac,chordae] = unpack_params(params); 

% total internal points in triangular domain 
total_internal = 3*N*(N+1)/2; 

% check if we have the with chordae version 
if isfield(params, 'chordae') && ~isempty(chordae)
    [m N_chordae] = size(chordae.C_left); 
else
    N_chordae = 0; 
end 

total_points = total_internal + 3*2*N_chordae; 

J = zeros(total_points,total_points); 

% always 6 pressure neighbors, which may or may not be in bounds
% relative indices of pressure here 
% numbered counter clockwise 
% ignore out of bounds indices, they are not in the pressure 
% bearing part of the surface 
pressure_nbrs = [ 0, -1; 
                  1, -1; 
                  1,  0; 
                  0,  1; 
                 -1,  1; 
                 -1,  0]';   


for j=1:N
    for k=1:N

        % in the triangle?
        if (j+k) < (N+2)

            % vertical offset does not change while differentiating this equation 
            range_current = linear_index_offset(j,k,N) + (1:3); 


            % pressure portion 
            % zero indexed loop because we are computing indices with mod n 
            for n=0:5

                j_nbr      = j + pressure_nbrs(1,mod(n  ,6)+1); 
                k_nbr      = k + pressure_nbrs(2,mod(n  ,6)+1); 
                j_nbr_next = j + pressure_nbrs(1,mod(n+1,6)+1); 
                k_nbr_next = k + pressure_nbrs(2,mod(n+1,6)+1);

                % if any index is zero, then 
                % the pressure term does not include this triangle
                if j_nbr_next && k_nbr_next && j_nbr && k_nbr

                    % Current has two terms from a product rule 
                    J(range_current, range_current) = J(range_current, range_current) ... 
                                                        - (p_0/6) * cross_matrix(X(:,j_nbr     ,k_nbr     ) - X(:,j,k)) ... 
                                                        + (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k)); 

                    % nbr term
                    % nbr gets differentiated away, and nbr_next stays 
                    % only added if this is internal 
                    if is_internal(j_nbr,k_nbr,N)
                        range_nbr       = linear_index_offset(j_nbr,k_nbr,N) + (1:3);
                        J(range_current, range_nbr) = J(range_current, range_nbr) - (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k));
                    end 

                    % nbr_next term
                    % nbr_next gets differentiated away, and nbr stays and gets a sign 
                    % only added if this is internal 
                    if is_internal(j_nbr_next,k_nbr_next,N)
                        range_nbr_next  = linear_index_offset(j_nbr_next,k_nbr_next,N) + (1:3);
                        J(range_current, range_nbr_next) = J(range_current, range_nbr_next) + (p_0/6) * cross_matrix(X(:,j_nbr,k_nbr) - X(:,j,k));    
                    end 

                end
            end 


            % u tension terms 
            for j_nbr = [j-1,j+1]

                k_nbr = k; 

                [X_nbr R_nbr idx_chordae left_side] = get_neighbor(params, filter_params, j_nbr, k_nbr); 

                J_tension = tension_jacobian(X(:,j,k),X_nbr,R(:,j,k),R_nbr,alpha,ref_frac); 

                % current term is always added in 
                % this gets no sign 
                % this is always at the current,current block in the matrix 
                J(range_current, range_current) = J(range_current, range_current) + J_tension; 

                % If the neighbor is an internal point, it also gets a Jacobian contribution 
                % This takes a sign
                if is_internal(j_nbr,k_nbr,N)
                    range_nbr  = linear_index_offset(j_nbr,k_nbr,N) + (1:3);
                    J(range_current, range_nbr) = J(range_current, range_nbr) - J_tension; 

                % If neighbor is on the chordae, it has a non zero index 
                % This is included here 
                elseif idx_chordae ~= 0
                    range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side); 
                    J(range_current, range_nbr) = J(range_current, range_nbr) - J_tension; 
                end 

            end 


            % v tension terms 
            for k_nbr = [k-1,k+1]

                j_nbr = j; 

                [X_nbr R_nbr idx_chordae left_side] = get_neighbor(params, filter_params, j_nbr, k_nbr); 

                J_tension = tension_jacobian(X(:,j,k),X_nbr,R(:,j,k),R_nbr,beta,ref_frac); 

                % current term is always added in 
                % this gets no sign 
                % this is always at the current,current block in the matrix 
                J(range_current, range_current) = J(range_current, range_current) + J_tension; 

                % If the neighbor is an internal point, it also gets a Jacobian contribution 
                % This takes a sign
                if is_internal(j_nbr,k_nbr,N)
                    range_nbr  = linear_index_offset(j_nbr,k_nbr,N) + (1:3);
                    J(range_current, range_nbr) = J(range_current, range_nbr) - J_tension; 

                % If neighbor is on the chordae, it has a non zero index 
                % This is included here 
                elseif idx_chordae ~= 0
                    range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side); 
                    J(range_current, range_nbr) = J(range_current, range_nbr) - J_tension; 
                end 

            end 

        end
    end
end




% chordae internal terms 
if N_chordae > 0
    
    [C_left, C_right, left_papillary, right_papillary, Ref_l, Ref_r, k_l, k_r, k_0, k_multiplier] = unpack_chordae(params.chordae); 
   
    for left_side = [true false];  
        
        if left_side
            C = C_left; 
            Ref = Ref_l; 
        else 
            C = C_right; 
            Ref = Ref_r; 
        end 
        
    
        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            % this is the same, updating the equations for this component 
            range_current = range_chordae(total_internal, N_chordae, i, left_side); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val j_nbr k_nbr] = get_nbr_chordae(params, i, nbr_idx, left_side); 

                % if the neighbor is in the chordae 
                if isempty(j_nbr) && isempty(k_nbr) 
                    range_nbr = range_chordae(total_internal, N_chordae, nbr_idx, left_side); 
                else
                    % neighbor is on the leaflet 
                    range_nbr = linear_index_offset(j_nbr,k_nbr,N) + (1:3);
                end 

                % tension Jacobian for this spring 
                J_tension = tension_jacobian(C(:,i), nbr, Ref(:,i), R_nbr, k_val, ref_frac); 

                % current always gets a contribution from this spring 
                J(range_current, range_current) = J(range_current, range_current) + J_tension; 

                % range may be empty if papillary muscle, in which case do nothing 
                if ~isempty(range_nbr)
                    J(range_current, range_nbr) = J(range_current, range_nbr) - J_tension; 
                end 
            end 
        end 
        
    end 
    
end 


nonzero_blocks = zeros(total_points/3); 

j = 1; 
for j_block=1:3:total_points
    k=1; 
    for k_block=1:3:total_points
        if J(j_block + (0:2), k_block + (0:2)) ~= zeros(3,3)
            nonzero_blocks(j,k) = 1; 
        end 
        k=k+1; 
    end
    j=j+1; 
end 

figure; 
spy(nonzero_blocks); 
title('single dot places at all nonzero blocks')
'whoa'











