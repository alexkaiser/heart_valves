function [leaflet valve] = set_tension_coeffs(leaflet, valve, tension_coeffs)
% 
% Returns three arrays with information about the geometry 
% 
% Output: 
%     is_internal          Boolean, true if 
%     is_bc                Point is a boundary condition that is fixed 
%     linear_idx_offset    In Jacobian, linear_idx_offset(j,k) + 1:3
%                          are the indices for the vector X(:,j,k)
% 


j_max                     = leaflet.j_max; 
k_min                     = leaflet.k_min; 
k_max                     = leaflet.k_max; 
n_rings_periodic          = leaflet.n_rings_periodic;
is_internal               = leaflet.is_internal; 

N_anterior                = valve.N_anterior; 
N_posterior               = valve.N_posterior;
commissural_leaflets      = valve.commissural_leaflets; 
if commissural_leaflets 
    N_commissure          = valve.N_commissure;  
else 
    N_commissure          = 0; 
end 

tension_base              = valve.p_physical / tension_coeffs.pressure_tension_ratio; 
dec_tension_coeff_base    = tension_coeffs.dec_tension_coeff_base; 
leaf_tension_base         = tension_coeffs.leaf_tension_base; 
root_tension_base         = tension_coeffs.root_tension_base ; 

alpha_anterior            = tension_coeffs.alpha_anterior  * tension_base; % circumferential 
beta_anterior             = tension_coeffs.beta_anterior   * tension_base; % radial
alpha_posterior           = tension_coeffs.alpha_posterior * tension_base; % circumferential 
beta_posterior            = tension_coeffs.beta_posterior  * tension_base; % radial
alpha_hoops               = tension_coeffs.alpha_hoops     * tension_base; % circumferential hoops  

c_circ_dec_anterior       = tension_coeffs.c_circ_dec_anterior       * dec_tension_coeff_base;  % circumferential 
c_rad_dec_anterior        = tension_coeffs.c_rad_dec_anterior        * dec_tension_coeff_base;  % radial
c_circ_dec_posterior      = tension_coeffs.c_circ_dec_posterior      * dec_tension_coeff_base;  % circumferential 
c_rad_dec_posterior       = tension_coeffs.c_rad_dec_posterior       * dec_tension_coeff_base;  % radial
c_circ_dec_hoops          = tension_coeffs.c_circ_dec_hoops          * dec_tension_coeff_base;  % radial hoops
c_rad_dec_hoops_anterior  = tension_coeffs.c_rad_dec_hoops_anterior  * dec_tension_coeff_base;  % radial hoops, anterior part 
c_rad_dec_hoops_posterior = tension_coeffs.c_rad_dec_hoops_posterior * dec_tension_coeff_base;  % radial hoops, posterior part 


if ~isfield(leaflet, 'chordae')
    error('Must initialize chordae prior to setting tensions'); 
end 


% coefficients for computing tensions 
alpha                     = zeros(j_max, k_max); 
beta                      = zeros(j_max, k_max);
c_dec_radial              = zeros(j_max, k_max); 
c_dec_circumferential     = zeros(j_max, k_max);


if commissural_leaflets 
    j_range_anterior   = (1:N_anterior); 
    j_range_right_comm = (1:N_commissure) + max(j_range_anterior); 
    j_range_posterior  = (1:N_posterior)  + max(j_range_right_comm); 
    j_range_left_comm  = (1:N_commissure) + max(j_range_posterior); 
    
    indices = [j_range_anterior, j_range_right_comm, j_range_posterior, j_range_left_comm]; 
    if ~all(indices == (1:j_max))
        error('Inconsistency in indices'); 
    end 
        
    alpha_commissure           = tension_coeffs.alpha_commissure           * tension_base;           % circumferential 
    beta_commissure            = tension_coeffs.beta_commissure            * tension_base;           % radial
    c_circ_dec_commissure      = tension_coeffs.c_circ_dec_commissure      * dec_tension_coeff_base; % circumferential 
    c_rad_dec_commissure       = tension_coeffs.c_rad_dec_commissure       * dec_tension_coeff_base; % radial 
    c_rad_dec_hoops_commissure = tension_coeffs.c_rad_dec_hoops_commissure * dec_tension_coeff_base; % radial in commissures 

else 
    j_range_anterior   = (1:N_anterior); 
    j_range_right_comm = [];
    j_range_posterior  = (1:N_posterior)  + max(j_range_anterior); 
    j_range_left_comm  = []; 
    
    indices = [j_range_anterior, j_range_posterior]; 
    if ~all(indices == (1:j_max))
        error('Inconsistency in indices'); 
    end 
    
end 


% circumferential hoops 
k_min_hoop = k_max - n_rings_periodic;

% radial anterior 
for j=j_range_anterior
    for k=k_min(j):(k_max-1)
        
        beta(j,k)         = beta_anterior; 
        
        if k < k_min_hoop
            c_dec_radial(j,k) = c_rad_dec_anterior; 
        else 
            c_dec_radial(j,k) = c_rad_dec_hoops_anterior; 
        end 
    end
end 

% radial posterior 
for j=j_range_posterior 
    for k=k_min(j):(k_max-1)
        
        beta(j,k) = beta_posterior; 

        if k < k_min_hoop
            c_dec_radial(j,k) = c_rad_dec_posterior; 
        else 
            c_dec_radial(j,k) = c_rad_dec_hoops_posterior; 
        end 
        
    end
end 

if commissural_leaflets 
    for j=j_range_right_comm
        for k=k_min(j):(k_max-1)
            
            beta(j,k)         = beta_commissure; 
        
            if k < k_min_hoop
                c_dec_radial(j,k) = c_rad_dec_commissure; 
            else 
                c_dec_radial(j,k) = c_rad_dec_hoops_commissure; 
            end 
        end
    end 

    % radial posterior 
    for j=j_range_left_comm 
        for k=k_min(j):(k_max-1)

            beta(j,k) = beta_commissure; 

            if k < k_min_hoop
                c_dec_radial(j,k) = c_rad_dec_commissure; 
            else 
                c_dec_radial(j,k) = c_rad_dec_hoops_commissure; 
            end 
        end
    end 
end 


% hoops in circumferential (hoop) direction 
for j=1:j_max 
    for k=k_min_hoop:(k_max-1)
        alpha(j,k) = alpha_hoops; 
        c_dec_circumferential(j,k) = c_circ_dec_hoops; 
    end 
end 

% cicumferential anterior 
% here need to ensure that neighbors are in bounds 
% also do not place a leaflet to leaflet circumferential spring
% only hoops connect leaflets 
for j=j_range_anterior(1:(end-1))
    
    % start at minimum, stop below hoop points 
    for k=k_min(j):(k_min_hoop-1)
        
        % spring is always owned by minimum neighbor 
        % j direction springs here 
        j_nbr = j + 1; 
        k_nbr = k; 
        
        if is_internal(j_nbr, k_nbr)
            alpha(j,k)                 = alpha_anterior; 
            c_dec_circumferential(j,k) = c_circ_dec_anterior; 
        end 
        
    end 
    
end 

% cicumferential posterior 
for j=j_range_posterior(1:(end-1))
    
    % start at minimum, stop below hoop points 
    for k=k_min(j):(k_min_hoop-1)
        
        % spring is always owned by minimum neighbor 
        % j direction springs here 
        j_nbr = j + 1; 
        k_nbr = k; 
        
        if is_internal(j_nbr, k_nbr)
            alpha(j,k) = alpha_posterior; 
            c_dec_circumferential(j,k) = c_circ_dec_posterior; 
        end 
        
    end 
    
end 

if commissural_leaflets 
    for j=j_range_right_comm(1:(end-1))

        % start at minimum, stop below hoop points 
        for k=k_min(j):(k_min_hoop-1)

            % spring is always owned by minimum neighbor 
            % j direction springs here 
            j_nbr = j + 1; 
            k_nbr = k; 

            if is_internal(j_nbr, k_nbr)
                alpha(j,k)                 = alpha_commissure; 
                c_dec_circumferential(j,k) = c_circ_dec_commissure; 
            end 

        end 

    end 

    % cicumferential posterior 
    for j=j_range_left_comm(1:(end-1))

        % start at minimum, stop below hoop points 
        for k=k_min(j):(k_min_hoop-1)

            % spring is always owned by minimum neighbor 
            % j direction springs here 
            j_nbr = j + 1; 
            k_nbr = k; 

            if is_internal(j_nbr, k_nbr)
                alpha(j,k) = alpha_commissure; 
                c_dec_circumferential(j,k) = c_circ_dec_commissure; 
            end 

        end 

    end 
end 


% set chordae tensions 
k_root                = tension_coeffs.k_root * tension_base * root_tension_base; 
k_0_1                 = tension_coeffs.k_0_1  * tension_base * leaf_tension_base;
chordae               = leaflet.chordae; 
c_dec_tension_chordae = tension_coeffs.c_dec_tension_chordae * dec_tension_coeff_base; 

for tree_idx = 1:leaflet.num_trees    

    free_edge_idx       = chordae(tree_idx).free_edge_idx; 
    
    n_leaves = size(free_edge_idx,1);  

    if (length(k_0_1) == 1) && (length(k_root) == 1)
        k_0_1_tmp   = k_0_1; 
        k_root_tmp  = k_root; 
    elseif (length(k_0_1) == leaflet.num_trees) && (length(k_root) == leaflet.num_trees)
        k_0_1_tmp   = k_0_1(tree_idx); 
        k_root_tmp  = k_root(tree_idx); 
    else
        error('Must supply values for all or values for exactly one')
    end 
    

    % Derived constants 
    
    % leaf force is total leaf force over number of leaves 
    k_0 = k_0_1_tmp / n_leaves; 
    
    % scaling formula on k_multiplier 
    % to achieve desired k_root 
    k_multiplier = 2.0 * (k_root_tmp/k_0_1_tmp)^(1/log2(n_leaves)); 
       
    % there are max_internal points on the tree 
    % leaves are not included as they are part of the leaflet 
    [m max_internal] = size(chordae(tree_idx).C); 
    
    % this parameter is the number of (not included) leaves in the tree
    NN = max_internal + 1; 
    
    % sanity check in building a balanced tree 
    n_tree = log2(NN);
    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 
    
    % Set tensions 
    
    % Each keeps parent wise spring constants 
    chordae(tree_idx).k_vals = zeros(max_internal,1); 

    % set spring constant data structures 
    num_at_level = n_leaves/2; 
    idx = max_internal; 

    % constants connecting to the leaflet are inherited
    % first internal constant to the tree is k_multiplier times that 
    k_running = k_multiplier*k_0; 
    while num_at_level >= 1

        for j=1:num_at_level
            chordae(tree_idx).k_vals(idx) = k_running; 
            idx = idx - 1; 
        end 

        k_running = k_running * k_multiplier; 
        num_at_level = num_at_level / 2; 
    end 

    % check that we actually got the desired root 
    tol = 1e-8; 
    if abs(chordae(tree_idx).k_vals(1) - k_root_tmp) > tol 
        error('Scaling incorrect at tree root, constants inconsistent'); 
    end 
        
    chordae(tree_idx).k_0                   = k_0; 
    chordae(tree_idx).k_multiplier          = k_multiplier;
    chordae(tree_idx).c_dec_tension_chordae = c_dec_tension_chordae(tree_idx); 
    
end 


% final scaling on target points 
valve.tension_base     = tension_base; 
valve.target_net       = valve.target_net_unscaled       * tension_base; 
valve.target_papillary = valve.target_papillary_unscaled * tension_base; 
valve.eta_net          = valve.eta_net_unscaled          * tension_base; 
valve.eta_papillary    = valve.eta_papillary_unscaled    * tension_base; 

% reset everything 
tension_coeffs.tension_base   = tension_base;
leaflet.chordae               = chordae;
leaflet.tension_coeffs        = tension_coeffs; 
leaflet.alpha                 = alpha; 
leaflet.beta                  = beta; 
leaflet.c_dec_radial          = c_dec_radial; 
leaflet.c_dec_circumferential = c_dec_circumferential; 
 
