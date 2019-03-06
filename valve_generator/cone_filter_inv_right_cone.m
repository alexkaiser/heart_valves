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

