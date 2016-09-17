function [F F_chordae_left F_chordae_right] = difference_equations(leaflet)
% 
% Evaluation of the global difference equations at j,k
% Requires reference configuration R 
% 
% Input
%     leaflet    Current parameters 
%
% Output
%     F        Values of all difference equation, 3 by triangular array 
% 

X                 = leaflet.X; 
R                 = leaflet.R; 
N                 = leaflet.N; 
p_0               = leaflet.p_0; 
alpha             = leaflet.alpha; 
beta              = leaflet.beta; 
ref_frac          = leaflet.ref_frac; 
C_left            = leaflet.chordae.C_left; 
C_right           = leaflet.chordae.C_right; 
Ref_l             = leaflet.chordae.Ref_l; 
Ref_r             = leaflet.chordae.Ref_r; 
chordae_idx_left  = leaflet.chordae_idx_left; 
chordae_idx_right = leaflet.chordae_idx_right;



[m N_chordae] = size(C_left); 

if leaflet.radial_and_circumferential
    error('Radial and circumferential fibers not implemented')
end 

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
        if leaflet.is_internal(j,k)

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
                
                if (j_nbr > 0) && (k_nbr > 0) && (leaflet.is_internal(j_nbr,k_nbr) || leaflet.is_bc(j_nbr,k_nbr))
                    
                    u_tangent_term = u_tangent_term + tension_linear(X(:,j,k),X(:,j_nbr,k_nbr),R(:,j,k),R(:,j_nbr,k_nbr),alpha,ref_frac) * (X(:,j_nbr,k_nbr) - X(:,j,k));
                
                end 
            end 
            
            v_tangent_term = zeros(3,1); 
            for k_nbr = [k-1,k+1]
                    
                j_nbr = j; 
                
                if (j_nbr > 0) && (k_nbr > 0) && (leaflet.is_internal(j_nbr,k_nbr) || leaflet.is_bc(j_nbr,k_nbr))
                
                    v_tangent_term = v_tangent_term + tension_linear(X(:,j,k),X(:,j_nbr,k_nbr),R(:,j,k),R(:,j_nbr,k_nbr),beta,ref_frac) * (X(:,j_nbr,k_nbr) - X(:,j,k));
                
                end 
            end 

            chordae_left_term = zeros(3,1); 
            if chordae_idx_left(j,k)
                    % current node has a chordae connection 
                    
                    kappa = leaflet.chordae.k_0; 
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet 
                    leaf_idx = chordae_idx_left(j,k) + N_chordae; 
                    
                    % then take the parent index of that number in chordae variables 
                    idx_chordae = floor(leaf_idx/2);  
                    
                    X_nbr = leaflet.chordae.C_left(:,idx_chordae); 
                    R_nbr = leaflet.chordae.Ref_l (:,idx_chordae);
                    
                    chordae_left_term = tension_linear(X(:,j,k),X_nbr,R(:,j,k),R_nbr,kappa,ref_frac) * (X_nbr - X(:,j,k));
                   
             end 
                
             chordae_right_term = zeros(3,1);    
             if chordae_idx_right(j,k)
                    % current node has a chordae connection 

                    kappa = leaflet.chordae.k_0; 
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet 
                    leaf_idx = chordae_idx_right(j,k) + N_chordae; 
                    
                    % then take the parent index of that number in chordae variables 
                    idx_chordae = floor(leaf_idx/2);  
                    
                    X_nbr = leaflet.chordae.C_right(:,idx_chordae); 
                    R_nbr = leaflet.chordae.Ref_r  (:,idx_chordae);
                    
                    chordae_right_term = tension_linear(X(:,j,k),X_nbr,R(:,j,k),R_nbr,kappa,ref_frac) * (X_nbr - X(:,j,k));

                end 
            
            
            F(:,j,k) = pressure_term + u_tangent_term + v_tangent_term + chordae_left_term + chordae_right_term;
            
        end
    end 
end 
    

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
            [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, nbr_idx, left_side); 

            tension = tension_linear(C(:,i), nbr, Ref(:,i), R_nbr, k_val, ref_frac) * (nbr - C(:,i));  

            if left_side
                F_chordae_left(:,i)  = F_chordae_left(:,i)  + tension; 
            else 
                F_chordae_right(:,i) = F_chordae_right(:,i) + tension; 
            end 

        end 

    end 
end 

