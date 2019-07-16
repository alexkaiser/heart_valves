function T_grad = tension_gradient_with_reference(X, X_nbr, R, k_spr, leaflet, collagen_constitutive)
% 
% Returns the tension in the linear constitutive law 
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

if exist('collagen_constitutive', 'var')
    collagen_constitutive_tmp = collagen_constitutive; 
elseif isfield(leaflet, 'collagen_constitutive') && leaflet.collagen_constitutive
    collagen_constitutive_tmp = leaflet.collagen_constitutive; 
else 
    collagen_constitutive_tmp = false; 
end 

if collagen_constitutive_tmp
    
    collagen_curve       = leaflet.collagen_curve; 
    a                    = collagen_curve.a; 
    b                    = collagen_curve.b; 
    full_recruitment     = collagen_curve.full_recruitment; 
    eta_collagen         = collagen_curve.eta_collagen; 
    
    E = norm(X - X_nbr)/R - 1.0; 
    
    if E < 0
        coeff = 0; 
    elseif E < full_recruitment
        coeff = k_spr * a * b * exp(b*E);
    else 
        coeff = k_spr * eta_collagen; 
    end 
    
    T_grad = -(coeff/R) * (X_nbr - X) / norm(X_nbr - X);
    
else 

    % default linear law 
    if ~isfield(leaflet, 'ref_frac')
        ref_frac = leaflet.ref_frac; 
    else 
        ref_frac = 1; 
    end 

    R      = ref_frac * R; 
    
    % Linear law has constant derivative 
    T_grad = -(k_spr/R) * (X_nbr - X) / norm(X_nbr - X);

end 
    
    