function J = tension_tangent_jacobian_with_reference(X, X_nbr, R, k, leaflet, collagen_constitutive)
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

if exist('collagen_constitutive', 'var')
    collagen_constitutive_tmp = collagen_constitutive; 
elseif isfield(leaflet, 'collagen_constitutive') && leaflet.collagen_constitutive
    collagen_constitutive_tmp = leaflet.collagen_constitutive; 
else 
    collagen_constitutive_tmp = false; 
end 

X_norm = norm(X_nbr - X);

% scalar tension gets differentiated
% gradient gets outer product with tangent 

tangent = (X_nbr - X)/X_norm; 
tension_gradient = tension_gradient_with_reference(X, X_nbr, R, k, leaflet, collagen_constitutive_tmp); 

% First term is outer product of gradient of tension with tangent 
J = tension_gradient * tangent'; 

% tangent gets differentiated 
J = J + tension_with_reference(X, X_nbr, R, k, leaflet, collagen_constitutive_tmp) * tangent_jacobian(X, X_nbr); 


