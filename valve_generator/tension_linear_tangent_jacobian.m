function J_tension = tension_linear_tangent_jacobian(X_current,X_nbr,R,k_spr,ref_frac)
% 
% Computes the contribution to the Jacobian for X and its neighbor
% 
% Input: 
%     X_current      Jacobian is taken with respect to this variable 
%     X_nbr          Relevant neighbor in X
%     R              Reference length 
%     k              Spring constant 
% 
% Output: 
%     J_tension      3x3 Jacobian for this tension term 
%                    Signs are NOT included, local will get a negative 
% 

    X_norm = norm(X_nbr - X_current);
    
    
    if ~exist('ref_frac', 'var')
        ref_frac = 1; 
    end 
    
    % jacobian is an outer product 
    % plus a multiple of the identity 
    J_tension = -k_spr / (X_norm^3) * (X_nbr - X_current)*((X_nbr - X_current)') ... 
                -tension_linear(X_current,X_nbr,R,k_spr,ref_frac)/X_norm * eye(3); 
     
end 

