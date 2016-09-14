function val = cone_filter_inv(X, leaflet)
%
% Inverts the cone filter at X
% 
% Input 
%     xi, eta 
%     leaflet   Struct with cone filter parameters 
%     
% Output
%     X               Vector of components 
% 

tol = 1e2 * eps;

a = leaflet.filter.a; 
r = leaflet.filter.r; 
h = leaflet.filter.h;

x = X(1); 
y = X(2); 
z = X(3); 

if z < 0.0
    error('cannot evaluate cone for negative z'); 
end 

% maximum angle in both coordinate systems 
theta_0 = acos(-a / sqrt(a^2 + r^2 + h^2));

% minimum angle in both coordinate systems to be on the front sheet
integrand = @(phi) - sqrt( (1 - (a/r) * sin(phi)).^2 + (h/r)^2) ./ (1 - (2*a/r) * sin(phi) + (a/r)^2 + (h/r)^2) ; 
theta_min = theta_0 + quadgk(integrand,0,pi/2,'RelTol',tol,'AbsTol',tol);


% rotate in to check if we are on the
psi = asin(r / sqrt(h^2 + r^2)); 
rot = [cos(psi), 0, -sin(psi); 0, 1, 0; sin(psi), 0, cos(psi)]; 
rotated = rot * X; 

% Are we anywhere near the triangle? 
% If we are on the triangle, the rotation should be 
if abs(rotated(1)) < tol 
    xi  = rotated(2); 
    eta = rotated(3); 
    


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
    
elseif (y >= 0) && (x >= 0) 
    
    % right cone, just use the inverse    
    val = cone_filter_inv_right_cone(X, leaflet); 
    
elseif (y >= 0) && (x < 0) 
    
    % reflect back to the right front and apply the map
    X(1) = -X(1); 
    val = cone_filter_inv_right_cone(X, leaflet); 
    
    % map the right front preimage to the right back preimage 
    val(1) = val(1) - a; 
    val = rotation_matrix(theta_min - pi/2) * val; 
    val(1) = -val(1); 
    val = rotation_matrix(-(theta_min - pi/2)) * val; 
    val(1) = val(1) + a; 
    
    
elseif (y < 0) && (x >= 0)     
    % left cone 
    
    % reflect to use the right cone inverse 
    X(2) = -X(2); 
    val = cone_filter_inv_right_cone(X, leaflet); 
    
    % and back to the left 
    val(1) = -val(1); 
    
elseif (y < 0) && (x < 0)     
    
    % reflect back to the left front 
    X(1) = -X(1);
    
    % now use the left front inverse
    X(2) = -X(2); 
    val = cone_filter_inv_right_cone(X, leaflet); 
    val(1) = -val(1); 
    
    % map the left front preimage to the leftå back preimage 
    val(1) = val(1) + a; 
    val = rotation_matrix(-(theta_min - pi/2)) * val; 
    val(1) = -val(1); 
    val = rotation_matrix(theta_min - pi/2) * val; 
    val(1) = val(1) - a; 
    
else 
    error('No valid parameter range found, out of range or bugs'); 
    
end 
    

































