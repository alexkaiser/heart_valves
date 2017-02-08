function [torus] = initialize_torus_data_structures(N, repulsive_potential)
% 
% Initializes data structures for full solve.  
% 
% Parameters are declared here.
% Should be a script, but want to return in the structures 
% 
% Input: 
%     N   Size parameter used throughout 
% 


% Main data structure with everything 
torus.N = N; 
torus.j_max = N; 
torus.k_max = N; 
torus.tol_global = 1e-10;
torus.max_it = 4000; 

torus.bead_slip = true; 
torus.leaflet_only = true; 

torus.repulsive_potential = false; 

if repulsive_potential
    error('repulsive potential not implemented for torus')
end 

torus.repulsive_power     = 1; 

% general mesh parameters 
torus.du = 1/N; 
torus.dv = 1/N; 

% fiber spacing util parameter 
dt = 1/N; 

% coefficient has units of 1/L, then gets squared 
c_0 = N;  
torus.repulsive_coeff  = 1e-7 * c_0^2; 

torus.diff_eqns = @difference_equations_torus; 
torus.jacobian  = @build_jacobian_torus; 
    
% major and minor radius of torus 
torus.r = 1; 
torus.R = 2;

% number of wraps in major and minor directions 

% minor direction may wrap zero times 
torus.n_wraps = 2; 

% major direction must wrap at least once 
torus.m_wraps = 1; 


% general solve parameters

% name 
torus.base_name = sprintf('torus_%d', N); 

% box width 
torus.L = 2.5; 

% pressure / spring constant ratio  
% ratio 6 is for N=32
% ratio = 6 seems to make everything very stiff 
% turn down by order of magnitude, see if it helps 
torus.pressure_tension_ratio = 1.5; 


% original spring constants were for N = 32 debug width
% spring constants get multiplied by 32/N, so they are halfed if N==64
% use this refintement number accordingly 
torus.refinement = N/32.0; 

torus.p_physical = 100; 

% scaling for target points 
torus.target_multiplier = 40; 

% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
% z direction is doubled 
torus.n_lagrangian_tracers = 8; 

% Uses configuration of X 
torus.X_config_is_reference = true; 

% places this many exact copies of the leaflet downward in z 
% spring constants are all reduced by num_copies 
% spacing is always half a mesh width 
torus.num_copies = 3; 

% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
torus.collagen_springs_leaflet = false; 


% Spring constants in two directions 
torus.alpha    =  1.0;  % circumferential 
torus.beta     =  1.0;  % radial 
torus.p_0      = -0.1;  % negative sign on anterior leaflet 


tor = @(u,v) [ sin(v) .* (torus.R + torus.r * cos(u)); ... 
               cos(v) .* (torus.R + torus.r * cos(u)); ...  
                                    torus.r * sin(u)]; 


torus.preimage = zeros(2,N,N);                                 
                                
torus.X = zeros(3,N,N);


for k=1:N
    
    % initial u conditions go from zero 
    u_0 = 2*pi*torus.r * (k-1) * dt; 
    
    for j=1:N
        
        u = 2*pi*torus.r*torus.n_wraps * dt * (j-1) + u_0; 
        v = 2*pi*torus.R*torus.m_wraps * dt * (j-1); 
                
        torus.preimage(:,j,k) = [u; v];  
        
        torus.X(:,j,k) = tor(u,v); 
    end 
end 

fig = figure; 
fig = torus_plot(torus, fig); 


'done with initialize'




fig = figure; 
hold on 

for k=1:N
    plot(torus.preimage(1,:,k), torus.preimage(2,:,k), 'o'); 
end 
axis equal 
title('preimage')


% set util arrays 
torus.is_internal       =  ones(torus.j_max, torus.k_max); 
torus.is_bc             = zeros(torus.j_max, torus.k_max); 
torus.linear_idx_offset = zeros(torus.j_max, torus.k_max); 
torus.point_idx_with_bc = zeros(torus.j_max, torus.k_max); 

count = 0; 
for k=1:torus.k_max
    for j=1:torus.j_max
        if torus.is_internal(j,k)
            torus.linear_idx_offset(j,k) = count; 
            count = count + 3; 
        end 
    end 
end

count = 0;
for k=1:torus.k_max
    for j=1:torus.j_max
        if torus.is_internal(j,k) || torus.is_bc(j,k)
            torus.point_idx_with_bc(j,k) = count; 
            count = count + 1; 
        end 
    end 
end











