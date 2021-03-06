function leaflet_with_reference = set_rest_lengths_and_constants(leaflet, valve)
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
chordae                = leaflet.chordae; 
chordae_idx            = leaflet.chordae_idx; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
du                     = leaflet.du; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
num_trees              = leaflet.num_trees; 
strain                 = valve.strain; 
% diastolic_increment    = valve.diastolic_increment; 

R_u = zeros(j_max, k_max); 
k_u = zeros(j_max, k_max); 
R_v = zeros(j_max, k_max); 
k_v = zeros(j_max, k_max); 


% allocation for tree leaflet connection points 
chordae_with_reference = chordae; 
for tree_idx = 1:num_trees 
    [m N_chordae] = size(chordae_with_reference(tree_idx).C);
    
    % zero spring constants 
    chordae_with_reference(tree_idx).k_vals = zeros(N_chordae,1); 
    
    % allocate rest lengths 
    chordae_with_reference(tree_idx).R_ch   = zeros(N_chordae,1); 
    
    % set parameters that should never be used to empty 
    chordae_with_reference(tree_idx).k_0 = []; 
    
    n_leaves = size(chordae_with_reference(tree_idx).free_edge_idx, 1); 
    
    if n_leaves ~= (N_chordae+1)
        error('Free edge and chordae sizes are inconsistent'); 
    end 
    
    chordae_with_reference(tree_idx).k_free_edge = zeros(n_leaves,1); 
    chordae_with_reference(tree_idx).R_free_edge = zeros(n_leaves,1); 
    
end 


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
        
        
        % current node has a chordae connection
        % data added to free edge array
        % 
        if chordae_idx(j,k).tree_idx
            if ~(is_internal(j,k) || is_bc(j,k))
                error('Chordae attachment to point that is neither internal nor boundary condition'); 
            end 
            
            tree_idx = chordae_idx(j,k).tree_idx; 

            [m N_chordae] = size(chordae(tree_idx).C);
            c_dec_tension_chordae = chordae(tree_idx).c_dec_chordae_leaf; 
            du_chordae = 1; 

            kappa = chordae(tree_idx).k_0;

            % index in current free edge array 
            i = chordae_idx(j,k).leaf_idx; 

            % index that free edge would have if on tree
            % remember that leaves are only in the leaflet
            leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

            % then take the parent index of that number in chordae variables
            idx_chordae = floor(leaf_idx/2);

            X_nbr = chordae(tree_idx).C(:,idx_chordae);
            tension = kappa; 

            if decreasing_tension && (kappa ~= 0)
                tension = tension + kappa * tension_decreasing(X, X_nbr, du_chordae, c_dec_tension_chordae) ; 
            end

            [chordae_with_reference(tree_idx).k_free_edge(i), ... 
             chordae_with_reference(tree_idx).R_free_edge(i)] ...
                = get_rest_len_and_spring_constants(X, X_nbr, tension, strain, leaflet); 

        end 
            
    end 
end


% chordae internal terms 
for tree_idx = 1:num_trees

    C = chordae(tree_idx).C; 
    [m N_chordae] = size(C);
    % c_dec_tension_chordae = chordae(tree_idx).c_dec_tension_chordae; 
        
    % normalize this, no mesh parameters in chordae computations 
    du_chordae = 1; 

    for i=1:N_chordae

        % only go parent wise here
        % child owns parent wise spring constants 
        parent = floor(i/2); 

        nbr_idx = parent;  
        
        % get the neighbors coordinates, reference coordinate and spring constants
        [nbr R_nbr k_val j_nbr k_nbr c_dec_tension_chordae] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 

        tension = k_val; 

        if decreasing_tension && (k_val ~= 0.0)
            tension = tension + k_val * tension_decreasing(C(:,i), nbr, du_chordae, c_dec_tension_chordae) ; 
        end

        [chordae_with_reference(tree_idx).k_vals(i), ...
         chordae_with_reference(tree_idx).R_ch(i)]  ... 
             = get_rest_len_and_spring_constants(C(:,i), nbr, tension, strain, leaflet); 
         
        % reset all papillary points according to diastolic skeleton 
        chordae_with_reference(tree_idx).root = chordae(tree_idx).root; %+ diastolic_increment; 
    end 
end 

% Copy all basic data structures 
leaflet_with_reference = leaflet; 

% Add new information 
leaflet_with_reference.R_u = R_u;
leaflet_with_reference.k_u = k_u;
leaflet_with_reference.R_v = R_v;
leaflet_with_reference.k_v = k_v;

leaflet_with_reference.chordae = chordae_with_reference; 

leaflet_with_reference.diff_eqns = @difference_equations_with_reference; 
leaflet_with_reference.jacobian  = @build_jacobian_with_reference;

leaflet_with_reference.ref_frac = 1.0; 

% % annular dilation if requested 
% if isfield(valve, 'dilation_radius') && (valve.dilation > 0) 
%     
%     r = valve.r; 
%     r_new = r + valve.dilation; 
%     
%     
% end 

if isfield(valve.skeleton, 'valve_ring_pts_early_systole') && isfield(valve.skeleton, 'papillary_early_systole')

    N_anterior       = valve.N_anterior; 
    N_posterior      = valve.N_posterior;
           
    % centered around zero 
    min_anterior = -leaflet.total_angle_anterior/2; 
    max_anterior =  leaflet.total_angle_anterior/2; 

    % mesh anterior inclusive of ends 
    mesh_anterior   = linspace(min_anterior, max_anterior, N_anterior); 

    % mesh posterior includes two anterior points
    min_anterior_wrapped = min_anterior + 2*pi; 
    mesh_posterior  = linspace(max_anterior, min_anterior_wrapped, N_posterior + 2);
    mesh_posterior  = mesh_posterior(2:(end-1)); 

    mesh = [mesh_anterior mesh_posterior];
    
    % reset bc values at ring 
    for j=1:j_max
        leaflet_with_reference.X(:,j,k_max) = interpolate_valve_ring_points(valve, mesh(j), valve.skeleton.valve_ring_pts_early_systole); 
    end 
            
    % reset papillary tip locations 
    for tree_idx = 1:leaflet_with_reference.num_trees
        % leaflet = add_chordae(leaflet, tree_idx); 
        leaflet_with_reference.chordae(tree_idx).root = valve.skeleton.papillary_early_systole(:,tree_idx); 
    end     
end 
