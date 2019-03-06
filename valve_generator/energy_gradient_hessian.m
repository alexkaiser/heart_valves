function [E, F, J] = energy_gradient_hessian(X_linearized, leaflet)
% 
% Returns energy, gradient and hessian for current leaflet data structures 
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

if ~isfield(leaflet, 'energy')
    error('Cannot run optimization without energy'); 
end

leaflet_copy = internal_points_to_2d(X_linearized, leaflet);

E = leaflet.energy(leaflet_copy); 

% Difference equations are gradient of energy 
if nargout > 1
    [F_2d F_chordae_left F_chordae_right] = leaflet.diff_eqns(leaflet_copy); 
    
    % Note that force is NEGATIVE gradient, but optimization just takes the gradient 
    F = -linearize_internal_points(leaflet, F_2d, F_chordae_left, F_chordae_right);
end

% Jacobian of gradient, Hessian of energy 
% Similarly, this is the second deriv of energy, so is NEGATIVE Jacobian of force
if nargout > 2
    J = -leaflet.jacobian(leaflet_copy); 
end

