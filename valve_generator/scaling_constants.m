
% Script of hacks on estimating good parameters with bad scaling 
% in form of good scaling. 
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


N = 32; 
n = log2(32); 

du_old = 1.274556735895149e-01; 
dv_old = 1.274556735895149e-01; 

alpha_old = 1; 
beta_old = 1; 

k_0_old = 1.8; 

k_multiplier_old = 1.8; 


k_root_old = 1.889568000000001e+01; 

% Since alpha, beta normalized, 
% alpha, beta, k_0 are 1,1,1.8 


% used directly with this number 
repulsive_coeff_old = 2.238985441466275e-03; 





% New version 

du_new = 1/32;  
dv_new = 1/32; 

alpha_new = 1; 
beta_new  = 1; 

k_m_new = 1.8; 

% leaflet tension coeffs 
% alpha/du 

% connections to free edge on LEAFLET part given 
% alpha * dv
% beta * du 

% same relationship says that k_0 should be mean, multiplied by 1.8 
% and scaled by du and dv 

k_0_32_new = 1.8 * (0.5) * (alpha_new * dv_new + beta_new * du_new)

% leaf base constant
% also equal to the total leaf force in the tree 
k_0_1_new = k_0_32_new * N 

% we know the k multiplier here, can use it to compute root tension 
k_root_new_from_multiplicative_scaling = (k_m_new)^n * k_0_32_new 

% generally k_m is computed to preserve total tension in the root and leaflet 
% apply the formula as a consistency check 
k_m_from_rule = 2.0 * (k_root_new_from_multiplicative_scaling / k_0_1_new)^(1/n)

% repulsive coeff used as 
% alpha/du * coeff * (du^2)
% repulsive_coeff_old = alpha_new/du_new * repulsive_coeff_new * (du_new^2)


repulsive_coeff_new = (du_new/alpha_new) * (1/du_new^2) * repulsive_coeff_old











