function J_tension = tension_linear_tangent_jacobian(X_current,X_nbr,R_current,R_nbr,k_spr,ref_frac)
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

    X_norm = norm(X_nbr - X_current);
    
    % jacobian is an outer product 
    % plus a multiple of the identity 
    J_tension = -k_spr / (X_norm^3) * (X_nbr - X_current)*((X_nbr - X_current)') ... 
                -tension_linear_over_norm(X_current,X_nbr,R_current,R_nbr,k_spr,ref_frac) * eye(3); 
     
end 

