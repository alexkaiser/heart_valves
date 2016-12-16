function J = replusive_jacobian(X_current,X_nbr,p)
% 
% Computes the local Jacobian block for the repulsive term
% 
% Input: 
%     X_current      Jacobian is taken with respect to this variable 
%     X_nbr          Relevant neighbor in X
%     p              Power in energy 
% 
% Output: 
%     J_tension      3x3 Jacobian for this location 
% 

    X_norm = norm(X_nbr - X_current);
    
    % jacobian is an outer product 
    % plus a multiple of the identity 
    J = -(p*(p+2)) / (X_norm^(p+4)) * (X_nbr - X_current)*((X_nbr - X_current)') + p / X_norm^(p+2) * eye(3); 
     
end 

