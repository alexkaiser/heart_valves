function val = cone_filter_inv(X, filter_params)
%
% Inverts the cone filter at X
% 
% Input 
%     xi, eta 
%     filter_params   Struct with cone filter parameters 
%     
% Output
%     X               Vector of components 
% 

tol = 1e2 * eps;

a = filter_params.a; 
r = filter_params.r; 
h = filter_params.h;

x = X(1); 
y = X(2); 
z = X(3); 

if z < 0.0
    error('cannot evaluate cone for negative z'); 
end 

if x < 0
    error('other sheet not yet implemented')
end 

% rotate in to check if we are on the
psi = asin(r / sqrt(h^2 + r^2)); 

rot = [cos(psi), 0, -sin(psi); 0, 1, 0; sin(psi), 0, cos(psi)]; 

rotated = rot * X; 

% Are we anywhere near the triangle? 
% If we are on the triangle, the rotation should be 
if abs(rotated(1)) < tol 
    xi  = rotated(2); 
    eta = rotated(3); 
    
    % maximum angle in both coordinate systems 
    theta_0 = acos(-a / sqrt(a^2 + r^2 + h^2)); 

    % only check if we are not on the bottom 
    if eta > 0
        theta_right = atan2(eta, xi - a);
        if theta_right < theta_0
            error('we are supposed to be in the triangle, but are on the right cone')
        end 

        theta_left = atan2(eta, -xi - a);
        if theta_left < theta_0
            error('we are supposed to be in the triangle, but are on the left cone')
        end    
    end 

    val = [xi; eta]; 
    
elseif y >= 0
    % right cone 
    val = cone_filter_inv_right_cone(X, filter_params); 
    
else 
    % left cone 
    
    % reflect to use the right cone inverse 
    X(2) = -X(2); 
    val = cone_filter_inv_right_cone(X, filter_params); 
    
    % and back to the left 
    val(1) = -val(1); 
end 
    





