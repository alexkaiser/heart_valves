function J = dec_tension_jacobian(X,X_nbr,du,c)
% 
% Computes the local Jacobian block for the decreasing tension term
% 
% Input: 
%     X_current      Jacobian is taken with respect to this variable 
%     X_nbr          Relevant neighbor in X
%     p              Power in energy 
% 
% Output: 
%     J_tension      3x3 Jacobian for this location 
% 

    X_norm = norm(X_nbr - X);
    
    % scalar tension gets differentiated
    % gradient gets outer product with tangent 
    
    outer_prod = (X_nbr - X)*((X_nbr - X)')/(X_norm^2); 
    scalar_derivative_term = - 2 * X_norm * (c * du)^(-2) / (1 +  (c * du)^(-2) * norm(X_nbr-X)^2)^2; 
     
    J = scalar_derivative_term * outer_prod; 
    
    % tangent gets differentiated 
    J = J + tension_decreasing(X,X_nbr,du,c) * tangent_jacobian(X,X_nbr); 
end 

