function [] = solve_valve_ibamr_output()



% note that if extra posterior is equal to commissure angle
% then the poster leaflet has angle pi 
% extra_posterior = pi/6; 

% each the posterior and anterior get this much extra
% overlap = pi/12; 


total_posterior = pi + pi/6 + pi/12; 

chordae_tree = true;

a = 1; 
r = 1.5; 
h = 3; 
N = 64; 

arbitrary_papillary_points = true; 
if arbitrary_papillary_points 
    % ct scan method 
    r = 1.606587877768772; 
    left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
    right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 
    h = 2; 
else 
    left_papillary  = [0; -a; 0]; 
    right_papillary = [0;  a; 0]; 
end 

% just send things more central to see what happens 
hack_anterior_papillary_placement = false; 
if hack_anterior_papillary_placement 
    left_papillary  = left_papillary  + [1; 0; 0]; 
    right_papillary = right_papillary + [1; 0; 0]; 
end 


min_angle_posterior = -total_posterior/2; 
max_angle_posterior =  total_posterior/2; 

filter_params_posterior.a = a; 
filter_params_posterior.r = r; 
filter_params_posterior.h = h;
filter_params_posterior.N = N;
filter_params_posterior.min_angle = min_angle_posterior;
filter_params_posterior.max_angle = max_angle_posterior;

% whole posterior leaflet is reflected after the computation
% papillary points must be reflected in x before the solve 
filter_params_posterior.left_papillary  = [-1; 1; 1] .* left_papillary; 
filter_params_posterior.right_papillary = [-1; 1; 1] .* right_papillary; 

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
    k_multiplier = 1.8;  % 2.0; 
    tree_frac = 0.5; 
    params_posterior = add_chordae(params_posterior, filter_params_posterior, k_0, k_multiplier, tree_frac, arbitrary_papillary_points); 
else 
    warning('Running without tree-based chordae currently untested a likely broken')
end 

fig = surf_plot(params_posterior, filter_params_posterior); 
title('Reference configuration of posterior surface'); 

'initial difference equation norm'
err = total_global_err(params_posterior, filter_params_posterior)

if chordae_tree 
    figure;
    'Jacobian build time = '
    tic; 
    J = build_jacobian(params_posterior, filter_params_posterior);
    toc; 
    spy(J); 
    title('initial nonzero pattern of jacobian')
end 


tol_global = 1e-10;

% when resolution increases 
% need to up the tolerance a little because of conditioning 
if N >= 256 
    tol_global = 1e-10;
end 

if N >= 512 
    tol_global = 1e-9;
end 

max_it_global = 40; 

plot_and_save_freq = 100; 
start_it = 0; 

err_over_time_posterior = zeros(max_it_global,1); 

if arbitrary_papillary_points
    % little higher to make a taller leaflet 
    h = 4; 
end 


% anterior 
total_anterior = pi; 
min_angle_anterior = -total_anterior/2; 
max_angle_anterior =  total_anterior/2;  

filter_params_anterior.a = a; 
filter_params_anterior.r = r; 
filter_params_anterior.h = h;
filter_params_anterior.N = N;
filter_params_anterior.min_angle = min_angle_anterior;
filter_params_anterior.max_angle = max_angle_anterior;

% anterior leaflet has no reflection
% adjust papillary coordinates slightly 
if arbitrary_papillary_points
    adjustment_anterior_x_coord = 0.0; 
    % This seems to be making problems, turn off for now
    % 0.1; 
else 
    adjustment_anterior_x_coord = 0.0; 
end 

filter_params_anterior.left_papillary  = left_papillary  + [adjustment_anterior_x_coord; 0; 0] ;
filter_params_anterior.right_papillary = right_papillary + [adjustment_anterior_x_coord; 0; 0] ;

% reference and initial surfaces are the same 
R_anterior = build_reference_surface(filter_params_anterior); 
X_anterior = R_anterior; 

params_anterior = pack_params(X_anterior,alpha,beta,N,p_0,R_anterior,ref_frac); 

if chordae_tree
    k_multiplier = 2.0; 
    params_anterior = add_chordae(params_anterior, filter_params_anterior, k_0, k_multiplier, tree_frac, arbitrary_papillary_points); 
end 

fig = surf_plot(params_anterior, filter_params_anterior); 
title('Reference configuration of anterior surface'); 

'initial difference equation norm'
err = total_global_err(params_anterior, filter_params_anterior)

err_over_time_anterior = zeros(max_it_global,1); 



% p_range = p_0; 
ref_frac_range = ref_frac; 

p_range = 0.0; %-(0:2.5:7.5); 


% if (N == 64) && arbitrary_papillary_points
%     ref_frac_range = [.6, ref_frac];  
% end 
% 
% if (N == 128) && (~arbitrary_papillary_points)
%     ref_frac_range = [.6,.65,.67,ref_frac];  
%     % with ref_frac_range = [.6,.65,.67,ref_frac];  
%     % this converges without any "scares" with highly increasing errors 
% end 
% 
% if (N == 128) && arbitrary_papillary_points
%     ref_frac_range = [.4, .5, .55, .6, .625, .65, .675, ref_frac];  
% end 
% 
% if (N == 256) && (~arbitrary_papillary_points)
%     ref_frac_range = [.5, .55, .6,.65,.67,ref_frac];  
% end 
% 
% if (N == 256) && arbitrary_papillary_points
%     ref_frac_range = [.3:.05:.55, .575, .6:.01:.67, .675:.005:.69, .692:.002:ref_frac];  
% end 


% debug values 
% p_range = -(0:2.5:5); 
% ref_frac_range = .5:.1:.6; 


for ref_frac = ref_frac_range
    
    % reset the reference fraction in updates 
    fprintf(1, 'Solving at ref_frac = %f\n\n' , ref_frac); 
    params_posterior.ref_frac = ref_frac; 
    params_anterior.ref_frac  = ref_frac; 
    
    for p_0 = p_range

        % when the pressure changes, just update the pressure and re-run the setup 
        % this is some type of continuation algorithm 
        params_posterior.p_0 = p_0; 
        params_anterior.p_0 = p_0; 

        [params_posterior pass err_over_time_posterior it] = solve_valve_auto_continuation(params_posterior, filter_params_posterior, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time_posterior, ref_frac, 'posterior'); 

        err_posterior = total_global_err(params_posterior, filter_params_posterior)

        % reflect posterior to actually have two leaflets
        params_posterior.X(1,:,:) = -params_posterior.X(1,:,:); 
        filter_params_posterior.left_papillary(1)  = -filter_params_posterior.left_papillary(1); 
        filter_params_posterior.right_papillary(1) = -filter_params_posterior.right_papillary(1); 
        
        if chordae_tree
            params_posterior.chordae.C_left(1,:,:)      = -params_posterior.chordae.C_left(1,:,:); 
            params_posterior.chordae.C_right(1,:,:)     = -params_posterior.chordae.C_right(1,:,:);
            params_posterior.chordae.left_papillary(1)  = -params_posterior.chordae.left_papillary(1);  
            params_posterior.chordae.right_papillary(1) = -params_posterior.chordae.right_papillary(1);  
        end 
        

        if pass 
            fprintf('Global solve passed posterior\n\n')
        else 
            fprintf('Global solve failed posterior\n\n')
            return; 
        end


        [params_anterior pass err_over_time_anterior it] = solve_valve_auto_continuation(params_anterior, filter_params_anterior, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time_anterior, ref_frac, 'anterior'); 

        
        err_anterior = total_global_err(params_anterior, filter_params_anterior)

        if pass 
            fprintf('Global solve passed anterior\n\n')
        else 
            fprintf('Global solve failed\n\n')
            return; 
        end        
        
        fig = figure; 
        fig = surf_plot(params_posterior, filter_params_posterior, fig);
        hold on 
        fig = surf_plot(params_anterior, filter_params_anterior, fig);
        title(sprintf('valve at p = %f', abs(p_0))); 

        % reflect back for solves 
        params_posterior.X(1,:,:) = -params_posterior.X(1,:,:); 
        filter_params_posterior.left_papillary(1)  = -filter_params_posterior.left_papillary(1); 
        filter_params_posterior.right_papillary(1) = -filter_params_posterior.right_papillary(1); 
        
        if chordae_tree
            params_posterior.chordae.C_left(1,:,:)      = -params_posterior.chordae.C_left(1,:,:); 
            params_posterior.chordae.C_right(1,:,:)     = -params_posterior.chordae.C_right(1,:,:);
            params_posterior.chordae.left_papillary(1)  = -params_posterior.chordae.left_papillary(1);  
            params_posterior.chordae.right_papillary(1) = -params_posterior.chordae.right_papillary(1);  
        end 
        
        
    end
end


% final reflection to ensure in right place 
params_posterior.X(1,:,:) = -params_posterior.X(1,:,:); 
filter_params_posterior.left_papillary(1)  = -filter_params_posterior.left_papillary(1); 
filter_params_posterior.right_papillary(1) = -filter_params_posterior.right_papillary(1); 

if chordae_tree
    params_posterior.chordae.C_left(1,:,:)      = -params_posterior.chordae.C_left(1,:,:); 
    params_posterior.chordae.C_right(1,:,:)     = -params_posterior.chordae.C_right(1,:,:);
    params_posterior.chordae.left_papillary(1)  = -params_posterior.chordae.left_papillary(1);  
    params_posterior.chordae.right_papillary(1) = -params_posterior.chordae.right_papillary(1);  
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

base_name = strcat(base_name, sprintf('_%d', N)); 

save(strcat(base_name, '_final_data')); 


rewrite_files_from_solved_valve(N); 












