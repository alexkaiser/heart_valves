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
    
    J = zeros(total_internal,total_internal); 
    
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
                    
                    [X_nbr R_nbr] = get_neighbor(params, filter_params, j_nbr, k_nbr); 
                    
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
                    end 
                    
                end 
                

                % v tension terms 
                for k_nbr = [k-1,k+1]
                    
                    j_nbr = j; 
                    
                    [X_nbr R_nbr] = get_neighbor(params, filter_params, j_nbr, k_nbr); 
                    
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
                    end 
                    
                end 
                
            end
        end
    end
end 


