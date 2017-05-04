function leaflet = get_util_arrays_bead_slip(leaflet, valve)
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
ring_k_idx                = leaflet.ring_k_idx; 
n_rings_periodic          = leaflet.n_rings_periodic;

N_anterior                = valve.N_anterior; 
N_posterior               = valve.N_posterior;
commissural_leaflets      = valve.commissural_leaflets; 
if commissural_leaflets 
    N_commissure          = valve.N_commissure;  
else 
    N_commissure          = 0; 
end 

tension_coeffs            = leaflet.tension_coeffs; 

alpha_anterior            = tension_coeffs.alpha_anterior;  % circumferential 
beta_anterior             = tension_coeffs.beta_anterior;   % radial
alpha_posterior           = tension_coeffs.alpha_posterior; % circumferential 
beta_posterior            = tension_coeffs.beta_posterior;  % radial
alpha_hoops               = tension_coeffs.alpha_hoops;     % circumferential hoops  

c_circ_dec_anterior       = tension_coeffs.c_circ_dec_anterior;  % circumferential 
c_rad_dec_anterior        = tension_coeffs.c_rad_dec_anterior   ;  % radial
c_circ_dec_posterior      = tension_coeffs.c_circ_dec_posterior ;  % circumferential 
c_rad_dec_posterior       = tension_coeffs.c_rad_dec_posterior  ;  % radial
c_circ_dec_hoops          = tension_coeffs.c_circ_dec_hoops     ;  % radial hoops
c_rad_dec_hoops_anterior  = tension_coeffs.c_rad_dec_hoops_anterior  ;  % radial hoops, anterior part 
c_rad_dec_hoops_posterior = tension_coeffs.c_rad_dec_hoops_posterior ;  % radial hoops, posterior part 

% data management stuff 
is_internal               = zeros(j_max, k_max); 
is_bc                     = zeros(j_max, k_max); 
linear_idx_offset         = zeros(j_max, k_max); 
point_idx_with_bc         = zeros(j_max, k_max); 

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



if leaflet.radial_and_circumferential 
    
    % radial and circumferential fiber layout 
    
    % valve ring at (j,ring_k_idx(j))
    for j=1:j_max 
        is_bc(j,ring_k_idx(j)) = true; 
    end 
    
    for j=1:j_max 
        
        % minimum k idx is variable 
        k = k_min(j);
        
        % keep going until hitting ring, which is always a bc 
        while ~is_bc(j,k)
            is_internal(j,k) = true; 
            k = k+1; 
        end 
    end 

else 
    error('diag fibers not implemented with bead slip')
end 


% free edge on leaflet only version is a b.c. point 
if leaflet.leaflet_only
    
    error('Not implemented in current version'); 
    
    for left_side = [true, false]
        if left_side
            free_edge_idx = leaflet.free_edge_idx_left; 
        else 
            free_edge_idx = leaflet.free_edge_idx_right; 
        end 

        for i=1:size(free_edge_idx, 1)
            j = free_edge_idx(i,1);
            k = free_edge_idx(i,2);
            is_internal(j,k) = false; 
            is_bc(j,k)       = true; 
        end        
    end 
end 

% Indices for Jacobian building 
count = 0; 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)
            linear_idx_offset(j,k) = count; 
            count = count + 3; 
        end 
    end 
end

% Indices for spring attachments 
count = 0;
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k) || is_bc(j,k)
            point_idx_with_bc(j,k) = count; 
            count = count + 1; 
        end 
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


for j=1:j_max 
    for k=k_min_hoop:(k_max-1)
        alpha(j,k) = alpha_hoops; 
        c_dec_circumferential(j,k) = c_circ_dec_hoops; 
    end 
end 

% radial anterior 
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

% radial posterior 
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


leaflet.is_internal           = is_internal;
leaflet.is_bc                 = is_bc;
leaflet.linear_idx_offset     = linear_idx_offset;
leaflet.point_idx_with_bc     = point_idx_with_bc;
leaflet.alpha                 = alpha; 
leaflet.beta                  = beta; 
leaflet.c_dec_radial          = c_dec_radial; 
leaflet.c_dec_circumferential = c_dec_circumferential; 

leaflet.j_range_anterior   = j_range_anterior; 
leaflet.j_range_right_comm = j_range_right_comm; 
leaflet.j_range_posterior  = j_range_posterior; 
leaflet.j_range_left_comm  = j_range_left_comm; 
