function Q = get_linear_tension_struct(X, X_nbr, R, R_nbr, kappa, ref_frac, j, k, j_nbr, k_nbr)
%
% Returns a struct for tension with the following fields
% 
%     val        Tension 
%     j          Current index 
%     k          Current index 
%     j_nbr      Neighboring index 
%     k_nbr      Neighboring index 
%     G          Jacobian, which is a gradient since tensions are scalar valued 
% 
% Input: 
%     X          Jacobian is taken with respect to this variable 
%     X_nbr      Relevant neighbor in X
%     R          Reference coordinate at current location 
%     R_nbr      Reference coordinate at neighbor location 
%     kappa      Spring constant 
%     ref_frac   Multiplier for rest length 
%     j          Current index 
%     k          Current index 
%     j_nbr      Neighboring index 
%     k_nbr      Neighboring index
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

X_norm =            norm(X - X_nbr);
R_norm = ref_frac * norm(R - R_nbr);

Q.val = kappa * (X_norm/R_norm - 1.0); 

Q.j = j; 
Q.k = k; 
Q.j_nbr = j_nbr; 
Q.k_nbr = k_nbr; 

Q.G = -(kappa/R_norm) * (X_nbr - X) / X_norm; 

