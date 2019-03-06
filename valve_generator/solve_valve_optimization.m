function [leaflet pass err E] = solve_valve_optimization(leaflet, tol)
% 
% Solves valve using built in matlab optimization routine 
% Leaflet must include an energy, difference equations and Jacobian 
% Note that difference equations are the gradient of energy
% and that the Jacobian is the Hessian of energy 
% 
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


energy_gradient_hessian_handle = @(X) energy_gradient_hessian(X, leaflet); 

% current configuration is initial guess 
X_linearized = linearize_internal_points(leaflet, leaflet.X, leaflet.chordae.C_left, leaflet.chordae.C_right); 

% should probably set some options here... 
options = optimset('LargeScale',        'on', ...
                   'Diagnostics',       'on', ...
                   'Display',           'iter-detailed', ...    % 'iter-detailed'
                   'GradObj',           'on', ...
                   'TolFun',             tol, ...
                   'TolX',               1.0e-12, ...
                   'Hessian',           'on'); 


[X_soln, E, exitflag, output, grad, hessian] = fminunc(energy_gradient_hessian_handle, X_linearized, options); 

'algorithm output:'
output
'algorithm output message:'
output.message


err = norm(grad); 

fprintf('flag = %d,\t E = %f,\t || F || = %f\n', exitflag, E, err); 


% send data back 
leaflet = internal_points_to_2d(X_soln, leaflet); 


% E_copy = leaflet.energy(leaflet); 
% err = total_global_err(leaflet); 
% 
% fprintf('E_copy = %f,\t || F || = %f\n', E_copy, err); 

pass = (err < tol); 









