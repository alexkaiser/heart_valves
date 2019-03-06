function X = cone_filter_right_cone(xi, eta, leaflet)
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

% Copyright (c) 2019, Alexander D. Kaiser
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

a = leaflet.filter.a; 
r = leaflet.filter.r; 
h = leaflet.filter.h; 

tol = 1e3 * eps;

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

if abs(theta - theta_0) > tol 

    % solve ODE for phi
    ode_solve_tol = 1e-10; 
    options = odeset('RelTol', ode_solve_tol, 'AbsTol', ode_solve_tol);
    
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

    if abs(z1 - z2) > 1e3*tol  
        error(sprintf('two formulas for z disagree, z1 = %f, z2 = %f', z1, z2)); 
    else 
        z = z1; 
    end 
end 

x = z*(r/h) * cos(phi); 
y = a*(1 - z/h) + z*(r/h) * sin(phi); 

X = [x;y;z]; 


