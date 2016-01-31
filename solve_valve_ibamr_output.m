function [] = solve_valve_ibamr_output()




% symmetric for now 
commissure_angle = 0; 

% note that if extra posterior is equal to commissure angle
% then the poster leaflet has angle pi 
extra_posterior = pi/6; 

% each the posterior and anterior get this much extra
overlap = pi/12; 

chordae_tree = true; 

a = 1; 
r = 1.5;
h = 3; 
N = 32; 

min_angle_posterior = -(pi/2 - commissure_angle/2 + extra_posterior/2 + overlap/2); 
max_angle_posterior =  (pi/2 - commissure_angle/2 + extra_posterior/2 + overlap/2); 

filter_params_posterior.a = a; 
filter_params_posterior.r = r; 
filter_params_posterior.h = h;
filter_params_posterior.N = N;
filter_params_posterior.min_angle = min_angle_posterior;
filter_params_posterior.max_angle = max_angle_posterior;

% reference and initial surfaces are the same 
R_posterior = build_reference_surface(filter_params_posterior); 
X_posterior = R_posterior; 

alpha     =  1.0; % spring constants in two directions 
beta      =  1.0;
p_0       = -0.0; 
ref_frac  =  0.7; 

params_posterior = pack_params(X_posterior,alpha,beta,N,p_0,R_posterior,ref_frac); 

if chordae_tree
    k_0 = 1; %0.000001; 
    k_multiplier = 1.75; 
    tree_frac = 0.5; 
    params_posterior = add_chordae(params_posterior, filter_params_posterior, k_0, k_multiplier, tree_frac); 
end 

fig = surf_plot(params_posterior, filter_params_posterior); 
title('Reference configuration of posterior surface'); 

'initial difference equation norm'
err = total_global_err(params_posterior, filter_params_posterior)

if chordae_tree 
    figure; 
    J = build_jacobian(params_posterior, filter_params_posterior); 
    spy(J); 
    title('initial nonzero pattern of jacobian')
end 


tol_global = 1e-11; 
max_it_global = 20; 

plot_and_save_freq = 10; 
start_it = 0; 

err_over_time_posterior = zeros(max_it_global,1); 


% anterior 
min_angle_anterior = -(pi/2 - commissure_angle/2 - extra_posterior/2 + overlap/2); 
max_angle_anterior =  (pi/2 - commissure_angle/2 - extra_posterior/2 + overlap/2); 

filter_params_anterior.a = a; 
filter_params_anterior.r = r; 
filter_params_anterior.h = h;
filter_params_anterior.N = N;
filter_params_anterior.min_angle = min_angle_anterior;
filter_params_anterior.max_angle = max_angle_anterior;

% reference and initial surfaces are the same 
R_anterior = build_reference_surface(filter_params_anterior); 
X_anterior = R_anterior; 

params_anterior = pack_params(X_anterior,alpha,beta,N,p_0,R_anterior,ref_frac); 

if chordae_tree
    params_anterior = add_chordae(params_anterior, filter_params_anterior, k_0, k_multiplier, tree_frac); 
end 

fig = surf_plot(params_anterior, filter_params_anterior); 
title('Reference configuration of anterior surface'); 

'initial difference equation norm'
err = total_global_err(params_anterior, filter_params_anterior)

err_over_time_anterior = zeros(max_it_global,1); 



% left commissural leaflet
% only fill in if commissural leaflet has nonzero angle 
if commissure_angle > 0.0
    % total hack on dimensions here for now     
    N_left = floor(N/2) + 1; 
    center = -(pi/2 - extra_posterior/2); 
    min_angle_left = center - commissure_angle/2; 
    max_angle_left = center + commissure_angle/2;
    left = true; 

    filter_params_left.a = a; 
    filter_params_left.r = r; 
    filter_params_left.h = h;
    filter_params_left.N = N_left;
    filter_params_left.min_angle = min_angle_left;
    filter_params_left.max_angle = max_angle_left;

    R_left = build_reference_surface_commissure(filter_params_left, left);
    X_left = R_left; 
    params_left = pack_params(X_left,alpha,beta,N_left,p_0,R_left,ref_frac); 

    fig = surf_plot_commissure(params_left, filter_params_left, left); 
    title('Reference configuration of surface'); 

    'initial difference equation norm'
    err_left = total_global_err_commissure(params_left, filter_params_left, left)

    err_over_time_left = zeros(max_it_global,1); 

    % right commissural leaflet
    % total hack on dimensions here for now 
    N_right = floor(N/2) + 1; 
    center = (pi/2 - extra_posterior/2); 
    min_angle_right = center - commissure_angle/2; 
    max_angle_right = center + commissure_angle/2;
    left = false; 

    filter_params_right.a = a; 
    filter_params_right.r = r; 
    filter_params_right.h = h;
    filter_params_right.N = N_right;
    filter_params_right.min_angle = min_angle_right;
    filter_params_right.max_angle = max_angle_right;

    R_right = build_reference_surface_commissure(filter_params_right, left);
    X_right = R_right; 
    params_right = pack_params(X_right,alpha,beta,N_right,p_0,R_right,ref_frac); 

    fig = surf_plot_commissure(params_right, filter_params_right, left); 
    title('Reference configuration of surface'); 

    'initial difference equation norm'
    err_right = total_global_err_commissure(params_right, filter_params_right, left)

    err_over_time_right = zeros(max_it_global,1); 
end 



% p_range = p_0; 
ref_frac_range = ref_frac; 

p_range = 0.0; %-(0:2.5:7.5); 
% ref_frac_range = .1:.1:1; 

% debug values 
% p_range = -(0:2.5:5); 
% ref_frac_range = .5:.1:.6; 


for ref_frac = ref_frac_range
    
    % reset the whole thing when the reference fraction changes
    params_posterior = pack_params(X_posterior,alpha,beta,N,p_0,R_posterior,ref_frac); 
    params_anterior = pack_params(X_anterior,alpha,beta,N,p_0,R_anterior,ref_frac); 
    
    if chordae_tree
        params_posterior = add_chordae(params_posterior, filter_params_posterior, k_0, k_multiplier, tree_frac); 
        params_anterior  = add_chordae(params_anterior,  filter_params_anterior , k_0, k_multiplier, tree_frac); 
    end 
    
    if commissure_angle > 0.0
        params_left = pack_params(X_left,alpha,beta,N_left,p_0,R_left,ref_frac); 
        params_right = pack_params(X_right,alpha,beta,N_right,p_0,R_right,ref_frac); 
        
        if chordae_tree
            error('chordae tree not yet supported for commissural leaflets'); 
        end 
        
    end 
    
    for p_0 = p_range

        % when the pressure changes, just update the pressure and re-run the setup 
        % this is some type of continuation algorithm 
        params_posterior.p_0 = p_0; 
        params_anterior.p_0 = p_0; 
        if commissure_angle > 0.0
            params_left.p_0 = p_0; 
            params_right.p_0 = p_0; 
        end
        



        [params_posterior pass err_over_time_posterior it] = solve_valve(params_posterior, filter_params_posterior, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time_posterior); 

        'difference equation norm'
        err_posterior = total_global_err(params_posterior, filter_params_posterior)

        % reflect posterior to actually have two leaflets
        params_posterior.X(1,:,:) = -params_posterior.X(1,:,:); 
        
        if chordae_tree
            params_posterior.chordae.C_left(1,:,:)  = -params_posterior.chordae.C_left(1,:,:); 
            params_posterior.chordae.C_right(1,:,:) = -params_posterior.chordae.C_right(1,:,:); 
        end 
        

        if pass 
            disp('Global solve passed posterior')
        else 
            disp('Global solve failed')
        end


        [params_anterior pass err_over_time_anterior it] = solve_valve(params_anterior, filter_params_anterior, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time_anterior); 

        'difference equation norm'
        err_anterior = total_global_err(params_anterior, filter_params_anterior)

        if pass 
            disp('Global solve passed anterior')
        else 
            disp('Global solve failed')
        end        
        
        if commissure_angle > 0.0
            left = true; 
            [params_left pass err_over_time_left it] = solve_commissure_leaflet(params_left, filter_params_left, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time_left, left); 


            if pass 
                disp('Global solve passed')
            else 
                disp('Global solve failed')
            end 

            'difference equation norm left'
            err_left = total_global_err_commissure(params_left, filter_params_left, left)


            left = false; 
            [params_right pass err_over_time_right it] = solve_commissure_leaflet(params_right, filter_params_right, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time_right, left); 


            if pass 
                disp('Global solve passed')
            else 
                disp('Global solve failed')
            end 

            'difference equation norm right'
            err_right = total_global_err_commissure(params_right, filter_params_right, left)
        end 

        
        fig = figure; 
        fig = surf_plot(params_posterior, filter_params_posterior, fig);
        hold on 
        fig = surf_plot(params_anterior, filter_params_anterior, fig);
        title(sprintf('valve at p = %f', abs(p_0))); 

        % reflect back for solves 
        params_posterior.X(1,:,:) = -params_posterior.X(1,:,:); 
        if chordae_tree
            params_posterior.chordae.C_left(1,:,:)  = -params_posterior.chordae.C_left(1,:,:); 
            params_posterior.chordae.C_right(1,:,:) = -params_posterior.chordae.C_right(1,:,:); 
        end
        
        
    end
end


% final reflection to ensure in right place 
params_posterior.X(1,:,:) = -params_posterior.X(1,:,:); 

if chordae_tree
    params_posterior.chordae.C_left(1,:,:)  = -params_posterior.chordae.C_left(1,:,:); 
    params_posterior.chordae.C_right(1,:,:) = -params_posterior.chordae.C_right(1,:,:); 
end


fig = figure; 
fig = surf_plot(params_posterior, filter_params_posterior, fig);
hold on 
fig = surf_plot(params_anterior, filter_params_anterior, fig);
title('final surface')




base_name = 'mitral'; 

if chordae_tree
    base_name = strcat(base_name, '_tree'); 
end 

L = 3; 
ratio = 6; 
p_physical = 100; 
target_multiplier = 100; 

% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
n_lagrangian_tracers = N / 2; 

output_to_ibamr_format(base_name, L, ratio, params_posterior, filter_params_posterior, params_anterior, p_physical, target_multiplier, n_lagrangian_tracers); 











