function J = tension_tangent_jacobian_with_reference(X, X_nbr, R, k, leaflet)
% 
% Computes the local Jacobian block for the tension term with reference config
% 
% Input: 
%     X_current      Jacobian is taken with respect to this variable 
%     X_nbr          Relevant neighbor in X
% 
% Output: 
%     J              3x3 Jacobian for this location 
% 

    X_norm = norm(X_nbr - X);
    
    % scalar tension gets differentiated
    % gradient gets outer product with tangent 
    
    tangent = (X_nbr - X)/X_norm; 
    tension_gradient = tension_gradient_with_reference(X, X_nbr, R, k, leaflet); 
     
    % First term is outer product of gradient of tension with tangent 
    J = tension_gradient * tangent'; 
    
    % tangent gets differentiated 
    J = J + tension_with_reference(X, X_nbr, R, k, leaflet) * tangent_jacobian(X, X_nbr); 
end 

