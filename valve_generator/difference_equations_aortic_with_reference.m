function F = difference_equations_aortic_with_reference(leaflet)
    % 
    % Evaluation of the global difference equations at j,k
    % Uses linear constitutive laws 
    % Requires rest lengths and spring constants to be previously
    % 
    % Input
    %     leaflet    Current parameters 
    %
    % Output
    %     F        Values of all difference equation, 3 by triangular array 
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

    X_current          = leaflet.X; 
    p_0                = leaflet.p_0; 
    j_max              = leaflet.j_max; 
    k_max              = leaflet.k_max; 
    is_internal        = leaflet.is_internal; 
    R_u                = leaflet.R_u;
    k_u                = leaflet.k_u;
    R_v                = leaflet.R_v;
    k_v                = leaflet.k_v;
    
    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
    end 
    
    if isfield(leaflet, 'targets_for_bcs') && leaflet.targets_for_bcs 
        targets_for_bcs = true; 
        k_target_net = leaflet.target_net; 
        k_target_papillary = leaflet.target_papillary; 
    else 
        targets_for_bcs = false; 
    end 
    
    if isfield(leaflet, 'target_length_check') && leaflet.target_length_check
        if isfield(leaflet, 'targets_for_bcs') && leaflet.targets_for_bcs
            if isfield(leaflet, 'ds')
                ds = leaflet.ds; 
                target_length_check = true; 
            else 
                warning('target_length_check is on, but ds is not, not checking')
                target_length_check = false; 
            end 
        else 
            warning('target_length_check is on, but targets_for_bcs is not, not checking'); 
            target_length_check = false; 
        end 
    else 
        target_length_check = false; 
    end 

    
    
    F_leaflet = zeros(size(X_current)); 

    % Internal leaflet part 
    for j=1:j_max
        for k=1:k_max
            if is_internal(j,k)

                X = X_current(:,j,k); 

                F_tmp = zeros(3,1);

                % pressure term first  
                if p_0 ~= 0
                    
                    [j_plus__1 j_minus_1 k_plus__1 k_minus_1 m] = get_pressure_nbrs(leaflet,j,k); 
                    
                    F_tmp = F_tmp + p_0 * m * cross(X_current(:,j_plus__1,k) - X_current(:,j_minus_1,k), X_current(:,j,k_plus__1) - X_current(:,j,k_minus_1));                     
                
                end 

                % u type fibers 
                for j_nbr_tmp = [j-1,j+1]
                    
                    k_nbr_tmp = k; 
                    
                    [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                    
                    if valid && (~target_spring) && (~target_k_no_j_spring)
                        X_nbr = X_current(:,j_nbr,k_nbr); 

                        tension = tension_with_reference(X, X_nbr, R_u(j_spr,k_spr), k_u(j_spr,k_spr), leaflet); 
                        F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 
                    
                    elseif valid && target_spring 
                        error('No j direction targets allowed'); 
                    end 
                    
                end 

                % v type fibers 
                for k_nbr_tmp = [k-1,k+1]

                    j_nbr_tmp = j; 
                    
                    [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                    
                    if valid && (~target_spring)
                        X_nbr = X_current(:,j_nbr,k_nbr); 
                        
                        tension = tension_with_reference(X, X_nbr, R_v(j_spr,k_spr), k_v(j_spr,k_spr), leaflet); 
                        F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 
                    
                    elseif valid && target_spring 
                        
                        if ~targets_for_bcs
                            error('Cannot ask for targets for bcs without flag set')
                        end
                        
                        X_nbr    = X_current(:,j_nbr,k_nbr);
                        tension_tangent = tension_zero_rest_length_linear_by_tangent(X, X_nbr, k_target_net); 
                        
                        if target_length_check
                            target_length = norm(X - X_nbr); 
                            if target_length > ds
                                fprintf('Found long target link in leaflet, L = %f, ds = %f\n', target_length, ds); 
                            end 
                        end
                        
                        % targets are absolute forces, no du here 
                        F_tmp = F_tmp + tension_tangent; 
                       
                    end 
                    
                end                
                
                F_leaflet(:,j,k) = F_tmp;

            end
        end 
    end
    
    F = linearize_internal_points(leaflet, F_leaflet); 
    
end 



