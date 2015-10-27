function X = cone_filter_right_cone(xi, eta, filter_params)
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
    X = [0;a;0]; 
    return; 
end 

% angle first 
theta = atan2(eta, xi - a);

dphi_dtheta = @(theta, phi) -(1 - (2*a/r) * sin(phi) + (a/r)^2 + (h/r)^2) ./ sqrt( (1 - (a/r) * sin(phi)).^2 + (h/r)^2); 

theta_0 = acos(-a / sqrt(a^2 + r^2 + h^2)); 

if (theta - theta_0) > tol 
    error('requesting theta greater than maximum allowed')
end 

integrand = @(phi) - sqrt( (1 - (a/r) * sin(phi)).^2 + (h/r)^2) ./ (1 - (2*a/r) * sin(phi) + (a/r)^2 + (h/r)^2) ; 
theta_min = theta_0 + quadgk(integrand,0,pi/2,'RelTol',tol,'AbsTol',tol); 
if (theta - theta_min) < -tol 
    error('requesting theta less than minimum allowed'); 
end 

theta_span = [theta_0, theta]; 
phi_0 = 0;  

if abs(theta - theta_0) > eps 

    % solve ODE for phi
    options = odeset('RelTol',tol,'AbsTol',tol);
    [theta_vec phi_vec] = ode45(dphi_dtheta, theta_span, phi_0, options); 
    phi = phi_vec(length(phi_vec)); 
else 
    % theta is the initial condition, so we know phi here without the ODE 
    phi = phi_0; 
end 

R = sqrt(r^2 - 2*a*r*sin(phi) + a^2 + h^2); 

% if norm([xi - a; eta]) > R 
%     error('attempting to evaluate at points which are off the sheet'); 
% end 

% catch for vertical points where first z formula does not work
if abs(theta - pi/2) < eps
    disp('hit vertical catch'); 
    z = h * eta / (R * sin(theta)); 
elseif (abs(theta - pi) < eps) || ((abs(theta) < eps))
    error('cannot take horizontal here'); 
else     
    z1 = h * (xi - a) / (R * cos(theta)); 
    z2 = h * eta / (R * sin(theta)); 

    if abs(z1 - z2) > 10*tol  
        error('two formulas for z disagree'); 
    else 
        z = z1; 
    end 
end 

x = z*(r/h) * cos(phi); 
y = a*(1 - z/h) + z*(r/h) * sin(phi); 

X = [x;y;z]; 


