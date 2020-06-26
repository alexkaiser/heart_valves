function [k R] = get_rest_len_and_spring_constants(X, X_nbr, tension, strain, leaflet, collagen_constitutive)
% 
% Given two points, a tension between them, and a desired specified strain 
% Compute rest length and spring constant that gives the current force 
% 
% Input: 
%     X          Point 
%     X_nbr      Other point 
%     tension    Force between the two (not that this must be an actual force, not a density)
%     strain     Desired strain 
% 
% Output: 
%     k          Spring constant k 
%                multiplies RELATIVE strain to get force 
%                This implies that k has units of force 
%     R          Rest length 
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

if strcmp(collagen_constitutive_tmp, 'aortic_circ')
    
    b = 57.456509400487398 * 2.0; 
    
    if strain >= 0
        % exponential through origin following May-Newman 2009 
        k = tension / (exp(b*strain) - 1.0); 
    else 
        error('Trying to set rest length with negative strain')
    end 
    
elseif strcmp(collagen_constitutive_tmp, 'aortic_rad')

    b = 22.397200094241359 * 2.0; 
        
    if strain >= 0
        % exponential through origin following May-Newman 2009 
        k = tension / (exp(b*strain) - 1.0); 
    else 
        error('Trying to set rest length with negative strain')
    end 

elseif strcmp(collagen_constitutive_tmp, 'linear_no_compressive')
 
    if strain >= 0
        k = tension / strain; 
    else 
        error('Trying to set rest length with negative strain')
    end  
    
elseif collagen_constitutive_tmp
    % mitral default 
    
    collagen_curve       = leaflet.collagen_curve; 
    a                    = collagen_curve.a; 
    b                    = collagen_curve.b; 
    full_recruitment     = collagen_curve.full_recruitment; 
    eta_collagen         = collagen_curve.eta_collagen;
    collagen_y_intercept = collagen_curve.collagen_y_intercept;
    
    if strain <= 0
        error('Cannot determine spring constant with negative strain.'); 
    
    elseif strain < full_recruitment
        
        % exponential part 
        k = tension / (a * (exp(b*strain) - 1)); 
    
    else
        
        % affine part of constitutive law
        k = tension / (eta_collagen * strain + collagen_y_intercept);
    
    end 
else 
    % linear law by default 
    k = tension / strain; 
end

% Rest length is determined by strain 
R = norm(X - X_nbr) / (strain + 1);

