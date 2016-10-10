function J = tangent_jacobian(X_current,X_nbr)
% 
% Computes the Jacobian of the unit tangent from X_current to X_nbr
% 
% Input: 
%     X_current      Jacobian is taken with respect to this variable 
%     X_nbr          Relevant neighbor in X
% 
% Output: 
%     J_tension      3x3 Jacobian for this unit tangent term 
% 

    X_norm = norm(X_nbr - X_current);
    
    % jacobian is an outer product 
    % plus a multiple of the identity 
    J = 1.0 / (X_norm^3) * (X_nbr - X_current)*((X_nbr - X_current)') - 1.0 / X_norm * eye(3); 
     
end 

