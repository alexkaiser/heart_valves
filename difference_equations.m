function F = difference_equations(params, filter_params)
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

left_papillary  = [0; -filter_params.a; 0]; 
right_papillary = [0;  filter_params.a; 0]; 

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
            

            % plus term always is included  
            u_tangent_term =  alpha * ( 1.0/(ref_frac*norm(R(:,j+1,k) - R(:,j,k))) - 1.0/norm(X(:,j+1,k) - X(:,j,k)) ) * (X(:,j+1,k) - X(:,j,k)); 
            
            % minus term may be a separate boundary condition 
            if j==1
                u_tangent_term = u_tangent_term - alpha * ( 1.0/(ref_frac*norm(R(:,j,k) - left_papillary)) - 1.0/norm(X(:,j,k) - left_papillary)) * (X(:,j,k) - left_papillary);
            else 
                u_tangent_term = u_tangent_term - alpha * ( 1.0/(ref_frac*norm(R(:,j,k) - R(:,j-1,k))) - 1.0/norm(X(:,j,k) - X(:,j-1,k)) ) * (X(:,j,k) - X(:,j-1,k)) ;
            end 
                
                
            % plus term always is included  
            v_tangent_term =   beta * ( 1.0/(ref_frac*norm(R(:,j,k+1) - R(:,j,k  ))) - 1.0/norm(X(:,j,k+1) - X(:,j,k  )) ) * (X(:,j,k+1) - X(:,j,k  )) ; 
            
            % minus term may be a separate boundary condition 
            if k==1
                v_tangent_term = v_tangent_term - beta * ( 1.0/(ref_frac*norm(R(:,j,k) - right_papillary)) - 1.0/norm(X(:,j,k) - right_papillary) ) * (X(:,j,k) - right_papillary);
            else 
                v_tangent_term = v_tangent_term - beta * ( 1.0/(ref_frac*norm(R(:,j,k) - R(:,j,k-1))) - 1.0/norm(X(:,j,k  ) - X(:,j,k-1)) ) * (X(:,j,k) - X(:,j,k-1));
            end 

            F(:,j,k) = pressure_term + u_tangent_term + v_tangent_term;
            
        end
    end 
end 
    


