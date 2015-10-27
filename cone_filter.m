function X = cone_filter(xi, eta, filter_params)
%
% Evaluates the cone filter at (xi, eta)
% 
% Input 
%     xi, eta         Euclidean coordinates in plane 
%     filter_params   Struct with cone filter parameters 
%     
% Output
%     X               Vector of components 
% 


a = filter_params.a; 
r = filter_params.r; 
h = filter_params.h; 

tol = 1e2 * eps;

% If on the pole return it
if abs(eta) < tol && abs(xi-a) < tol 
    X = [0; a;0]; 
    return; 
end 

% If on the other pole return it
if abs(eta) < tol && abs(xi+a) < tol 
    X = [0;-a;0]; 
    return; 
end 

% maximum angle in both coordinate systems 
theta_0 = acos(-a / sqrt(a^2 + r^2 + h^2)); 

% angle first 
theta_right = atan2(eta, xi - a);

% if angle is less than the leftmost ray 
% then we are on the right sheet 
if theta_right <= theta_0
    X = cone_filter_right_cone(xi, eta, filter_params); 
    return 
end 

% Other angle 
% use reflected coordinates, (-xi,eta)
% and the (a,0) polar coordinate origin 
theta_left = atan2(eta, -xi - a);

if theta_left <= theta_0
    % apply the reflection in xi to get the corresponding preimage 
    X = cone_filter_right_cone(-xi, eta, filter_params); 
    
    % reflect back in y in the 3d setting  
    X(2) = -X(2); 
    return; 
end 

% we are on the triangle, rotate out 
psi = asin(r / sqrt(h^2 + r^2)); 

X = zeros(3,1); 
X(1) = sin(psi) * eta; 
X(2) =             xi; 
X(3) = cos(psi) * eta; 

