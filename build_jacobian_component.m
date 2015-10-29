function J = build_jacobian_component(params, filter_params,j,k)
% 
% Builds the Jacobian for the current index and parameter values 
% 
% Input 
%      params    Current parameter values
%      j,k       Indices, must be internal to the arrays
% 
% Output 
%      J         Jacobian of difference equations 


    [X,alpha,beta,N,p_0,R] = unpack_params(params); 

    left_papillary  = [0; -filter_params.a; 0]; 
    right_papillary = [0;  filter_params.a; 0]; 

%    F = zeros(size(X)); 

    % total internal points in triangular domain 
    % total_internal = 3*N*(N+1)/2; 
    
    J = zeros(3,3); 

    
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
                

    % in the triangle?
    if (j+k) < (N+2)

        % vertical offset does not change while differentiating this equation 
        % range_current = linear_index_offset(j,k,N) + (1:3); 


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
                J = J - (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k)) ... 
                      + (p_0/6) * cross_matrix(X(:,j_nbr     ,k_nbr     ) - X(:,j,k)); 


% 
% No neighbors for these terms 
% 
                                                
%                 % nbr_next term
%                 % nbr_next gets differentiated away, and nbr stays and gets a sign 
%                 % only added if this is internal 
%                 if is_internal(j_nbr_next,k_nbr_next,N)
%                     range_nbr_next  = linear_index_offset(j_nbr_next,k_nbr_next,N) + (1:3);
%                     J(range_current, range_nbr_next) = J(range_current, range_nbr_next) - (p_0/6) * cross_matrix(X(:,j_nbr,k_nbr) - X(:,j,k));    
%                 end 
% 
% 
%                 % nbr term
%                 % nbr gets differentiated away, and nbr_next stays 
%                 % only added if this is internal 
%                 if is_internal(j_nbr,k_nbr,N)
%                     range_nbr       = linear_index_offset(j_nbr,k_nbr,N) + (1:3);
%                     J(range_current, range_nbr) = J(range_current, range_nbr) - (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k));
%                 end 

            end
        end 




        % u tension terms 
        for j_nbr = [j-1,j+1]

            k_nbr = k; 

            if j_nbr == 0
                X_nbr = left_papillary;
                R_nbr = left_papillary;
            else 
                X_nbr = X(:,j_nbr,k_nbr); 
                R_nbr = R(:,j_nbr,k_nbr); 
            end 

            J_tension = tension_term(X(:,j,k),X_nbr,R(:,j,k),R_nbr,alpha); 

            % current term is always added in 
            % this gets a negative sign 
            % this is always at the current,current block in the matrix 
            J = J - J_tension; 

%             % If the neighbor is an internal point, it also gets a Jacobian contribution 
%             % This takes a positive sign
%             if is_internal(j_nbr,k_nbr,N)
%                 range_nbr  = linear_index_offset(j_nbr,k_nbr,N) + (1:3);
%                 J(range_current, range_nbr) = J(range_current, range_nbr) + J_tension; 
%             end 

        end 


        % v tension terms 
        for k_nbr = [k-1,k+1]

            j_nbr = j; 

            if k_nbr == 0
                X_nbr = right_papillary;
                R_nbr = right_papillary;
            else 
                X_nbr = X(:,j_nbr,k_nbr); 
                R_nbr = R(:,j_nbr,k_nbr); 
            end 

            J_tension = tension_term(X(:,j,k),X_nbr,R(:,j,k),R_nbr,beta); 

            % current term is always added in 
            % this gets a negative sign 
            % this is always at the current,current block in the matrix 
            J = J - J_tension; 

%             % If the neighbor is an internal point, it also gets a Jacobian contribution 
%             % This takes a positive sign
%             if is_internal(j_nbr,k_nbr,N)
%                 range_nbr  = linear_index_offset(j_nbr,k_nbr,N) + (1:3);
%                 J(range_current, range_nbr) = J(range_current, range_nbr) + J_tension; 
%             end 

        end 

    end
    
end 


function val = is_internal(j,k,N)
%  
% Checks whether a given index is an internal point 
% 

    if (j < 1) || (k < 1)
        val = false; 
        return 
    elseif (j+k) >= (N+2)   
        val = false; 
        return; 
    end 
    
    val = true; 
end 


function idx = linear_index_offset(j,k,N)
    %
    % Maps j,k in 3d to correct offset in flattened array 
    % Assumse N internal points with triangle domain 
    % 
      
    prev_values = N*(N+1)/2 - (N-k+1)*(N-k+2)/2; 
    idx = 3 * (j-1 + prev_values); 
end 


function J_tension = tension_term(X_current,X_nbr,R_current,R_nbr,k_spr)
% 
% Computes the contribution to the Jacobian for X and its neighbor
% 
% Input: 
%     X_current      Jacobian is taken with respect to this variable 
%     X_nbr          Relevant neighbor in X
%     R_current      Reference coordinate at current location 
%     R_nbr          Reference coordinate at neighbor location 
%     k              Spring constant 
% 
% Output: 
%     J_tension      3x3 Jacobian for this tension term 
%                    Signs are NOT included, local will get a negative 
% 

    J_tension = zeros(3,3);

    % u difference norms which are used repeatedly
    X_norm = norm(X_nbr - X_current);
    R_norm = norm(R_nbr - R_current);

   
    for l=1:3
        for m=1:3

            J_tension(l,m) = (X_nbr(l) - X_current(l)) * (X_nbr(m) - X_current(m)) / (X_norm^3) ; 
                    
            % diagonal term has an extra term
            if l == m
                J_tension(l,m) = J_tension(l,m) + (1.0/R_norm - 1.0/X_norm);                       
            end
        end
    end
    
    J_tension = k_spr * J_tension;
    
end 

function C = cross_matrix(x)
%
% Returns the 3x3 matrix which has the action of applying x cross 
% 
    C = [   0   x(3) -x(2); 
         -x(3)    0   x(1); 
          x(2) -x(1)    0];
      
end 





