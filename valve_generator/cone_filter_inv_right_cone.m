function val = cone_filter_inv_right_cone(X, leaflet)
%
% Inverts the cone filter at X
% 
% Input 
%     xi, eta 
%     Leaflet    Struct with cone filter parameters 
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

if z <= 0.0
    error('cannot evaluate cone for negative z'); 
end 
    
phi_1 = acos(  x     *h/(z*r)); 
phi_2 = asin( (y - a)*h/(z*r) + a/r); 
if abs(phi_1 - phi_2) > tol;  
    phi_1
    phi_2
    diff = abs(phi_1 - phi_2)
    warning('two formulas for phi disagree'); 
    phi = phi_1; 
else 
    phi = phi_1;
end 

if (phi < 0) || (phi > pi/2)
    error('phi out of rance'); 
end 

R = sqrt(r^2 - 2*a*r*sin(phi) + a^2 + h^2); 

theta_0 = acos(-a / sqrt(a^2 + r^2 + h^2)); 

integrand = @(phi) - sqrt( (1 - (a/r) * sin(phi)).^2 + (h/r)^2) ./ (1 - (2*a/r) * sin(phi) + (a/r)^2 + (h/r)^2) ; 


theta = theta_0 + quadgk(integrand,0,phi,'RelTol',tol,'AbsTol',tol); 

xi  = a + (R*z/h)*cos(theta); 
eta =     (R*z/h)*sin(theta); 

val = [xi; eta]; 

