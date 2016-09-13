function [F F_chordae_left F_chordae_right] = difference_equations(params, filter_params)
% 
% Evaluation of the global difference equations at j,k
% Uses 
% Requires reference configuration R 
% 
% Input
%     params   Current parameters 
%     filter_params  Cone filter paramaters 
%
% Output
%     F        Values of all difference equation, 3 by triangular array 
% 


[X,alpha,beta,N,p_0,R,ref_frac] = unpack_params(params); 

F = zeros(size(X)); 

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

            % pressure term first  
            pressure_term = zeros(3,1); 
            
            if p_0 ~= 0
                % zero indexed loop because we are computing indices with mod n 
                for n=0:5

                    j_nbr      = j + pressure_nbrs(1,mod(n  ,6)+1);
                    k_nbr      = k + pressure_nbrs(2,mod(n  ,6)+1);
                    j_nbr_next = j + pressure_nbrs(1,mod(n+1,6)+1);
                    k_nbr_next = k + pressure_nbrs(2,mod(n+1,6)+1);

                    % if any index is zero, then 
                    % the pressure term does not include this value
                    if j_nbr_next && k_nbr_next && j_nbr && k_nbr
                        pressure_term = pressure_term + (p_0/6) * cross(X(:,j_nbr,k_nbr) - X(:,j,k), X(:,j_nbr_next,k_nbr_next) - X(:,j,k));                     
                    end 

                end 
            end 
            
            u_tangent_term = zeros(3,1); 
            for j_nbr = [j-1,j+1]
                    
                k_nbr = k; 

                [X_nbr R_nbr] = get_neighbor(params, filter_params, j_nbr, k_nbr); 
                
                u_tangent_term = u_tangent_term + tension_linear(X(:,j,k),X_nbr,R(:,j,k),R_nbr,alpha,ref_frac) * (X_nbr - X(:,j,k));

            end 
            
            v_tangent_term = zeros(3,1); 
            for k_nbr = [k-1,k+1]
                    
                j_nbr = j; 
                
                [X_nbr R_nbr] = get_neighbor(params, filter_params, j_nbr, k_nbr); 
                
                v_tangent_term = v_tangent_term + tension_linear(X(:,j,k),X_nbr,R(:,j,k),R_nbr,beta,ref_frac) * (X_nbr - X(:,j,k));
            end 

            F(:,j,k) = pressure_term + u_tangent_term + v_tangent_term;
            
        end
    end 
end 
    

% additional tension terms for chordae if appropriate 
if (~isfield(params, 'chordae')) || isempty(params.chordae)
    F_chordae_left  = 0; 
    F_chordae_right = 0; 
else 
    
    [C_left, C_right, left_papillary, right_papillary, Ref_l, Ref_r, k_l, k_r, k_0, k_multiplier] = unpack_chordae(params.chordae); 
    
    F_chordae_left  = zeros(size(C_left )); 
    F_chordae_right = zeros(size(C_right)); 
    
    [m N_chordae] = size(C_left); 
    
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

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val] = get_nbr_chordae(params, i, nbr_idx, left_side); 

                tension = tension_linear(C(:,i), nbr, Ref(:,i), R_nbr, k_val, ref_frac) * (nbr - C(:,i));  

                if left_side
                    F_chordae_left(:,i)  = F_chordae_left(:,i)  + tension; 
                else 
                    F_chordae_right(:,i) = F_chordae_right(:,i) + tension; 
                end 
                
            end 

        end 
    end 
        
end 
