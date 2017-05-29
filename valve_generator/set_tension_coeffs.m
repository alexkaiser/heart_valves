function leaflet = set_tension_coeffs(leaflet, valve, tension_coeffs)
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
        
    alpha_commissure           = tension_coeffs.alpha_commissure;      % circumferential 
    beta_commissure            = tension_coeffs.beta_commissure;       % radial
    c_circ_dec_commissure      = tension_coeffs.c_circ_dec_commissure; % circumferential 
    c_rad_dec_commissure       = tension_coeffs.c_rad_dec_commissure;  % radial 
    c_rad_dec_hoops_commissure = tension_coeffs.c_rad_dec_hoops_commissure; % radial in commissures 

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

leaflet.tension_coeffs        = tension_coeffs; 
leaflet.alpha                 = alpha; 
leaflet.beta                  = beta; 
leaflet.c_dec_radial          = c_dec_radial; 
leaflet.c_dec_circumferential = c_dec_circumferential; 
 
