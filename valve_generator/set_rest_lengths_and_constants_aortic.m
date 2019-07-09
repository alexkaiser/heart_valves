function leaflet_with_reference = set_rest_lengths_and_constants_aortic(leaflet, valve)
% 
% Assignes spring constants and rest lengths such that the current 
% valve configuration has uniform strain as specified here 
% 
% Tensions are computed from the configuration 
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

X_current              = leaflet.X; 
alpha                  = leaflet.alpha; 
beta                   = leaflet.beta; 
c_dec_radial           = leaflet.c_dec_radial; 
c_dec_circumferential  = leaflet.c_dec_circumferential; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
du                     = leaflet.du; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
strain                 = valve.strain; 
% diastolic_increment    = valve.diastolic_increment; 

R_u = zeros(j_max, k_max); 
k_u = zeros(j_max, k_max); 
R_v = zeros(j_max, k_max); 
k_v = zeros(j_max, k_max); 

if isfield(leaflet, 'decreasing_tension') && leaflet.decreasing_tension
    decreasing_tension = true; 
else 
    decreasing_tension = false; 
end 
    

% Internal leaflet part 
for j=1:j_max
    for k=1:k_max
        
        X = X_current(:,j,k); 
        
        if is_internal(j,k)    

            % u type fibers 
            for j_nbr_tmp = [j-1,j+1]

                k_nbr_tmp = k; 

                [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                if valid               
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    alpha_tmp     = alpha(j_spr,k_spr); 
                    c_dec_tension = c_dec_circumferential(j_spr,k_spr); 

                    tension = alpha_tmp;  

                    if decreasing_tension && (alpha_tmp ~= 0)
                        tension = tension + alpha_tmp * tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                    end 

                    % Here tension is a force per unit length 
                    % So we must multiply by a legnth element to get force 
                    % Take the opposing length element 
                    tension = du * tension; 

                    [k_u(j_spr,k_spr) R_u(j_spr,k_spr)] = get_rest_len_and_spring_constants(X, X_nbr, tension, strain, leaflet); 
                
                end 
                
            end 


            % v type fibers 
            for k_nbr_tmp = [k-1,k+1]

                j_nbr_tmp = j; 

                [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                if valid
                    X_nbr = X_current(:,j_nbr,k_nbr); 
                    
                    beta_tmp      = beta(j_spr,k_spr); 
                    c_dec_tension = c_dec_radial(j_spr,k_spr); 

                    tension = beta_tmp; 

                    if decreasing_tension && (beta_tmp ~= 0)
                        tension = tension + beta_tmp * tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                    end

                    tension = du * tension; 

                    [k_v(j_spr,k_spr) R_v(j_spr,k_spr)] = get_rest_len_and_spring_constants(X, X_nbr, tension, strain, leaflet); 
                
                end 
                
            end 
        end
            
    end 
end

% Copy all basic data structures 
leaflet_with_reference = leaflet; 

% Add new information 
leaflet_with_reference.R_u = R_u;
leaflet_with_reference.k_u = k_u;
leaflet_with_reference.R_v = R_v;
leaflet_with_reference.k_v = k_v;

leaflet_with_reference.diff_eqns = @difference_equations_aortic_with_reference; 
leaflet_with_reference.jacobian  = @build_jacobian_aortic_with_reference;

leaflet_with_reference.ref_frac = 1.0; 
