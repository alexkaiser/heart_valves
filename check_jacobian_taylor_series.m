% check_jacobian_taylor_series


% checks that the taylor series, using the jacobian included, really works 





epsilon_vals = 10.^(-1:-1:-8); 

errors = zeros(size(epsilon_vals)); 


a = 1; 
r = 1.5;
h = 2; 
N = 4; 
min_angle = -pi/2; 
max_angle =  pi/2; 

filter_params.a = a; 
filter_params.r = r; 
filter_params.h = h;
filter_params.N = N;
filter_params.min_angle = min_angle;
filter_params.max_angle = max_angle;

% reference and initial surfaces are the same 
R = build_reference_surface(filter_params); 

X = R; 
alpha     =  1.0; % spring constants in two directions 
beta      =  1.0;
p_0       = -2.0; 
ref_frac  =  0.5; 

params = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

chordae_tree = true; 

if chordae_tree
    k_0 = 1; 
    k_multiplier = 2; 
    tree_frac = 0.5; 
    params = add_chordae(params, filter_params, k_0, k_multiplier, tree_frac); 
    chordae = params.chordae;
    [m N_chordae] = size(params.chordae.C_left); 
else 
    chordae = []; 
    params = pack_params(X,alpha,beta,N,p_0,R,ref_frac,chordae); 
end 


% difference eqns at center do not change 
F = difference_equations(params, filter_params); 
F_linearized = linearize_internal_points(F, params); 

% jacobian does not change 
J = build_jacobian(params, filter_params); 

figure; 
spy(J); 
title('jacobian nonzero structure in jacobian tester')

% perturbation also does not change 
Z = zeros(size(X)); 
for j=1:params.N
    for k=1:params.N
        % in the triangle?
        if (j+k) < (params.N+2)
            Z(:,j,k) = rand(3,1);  
        end 
    end 
end 

params_Z = pack_params(Z,alpha,beta,N,p_0,R,ref_frac,chordae); 

if chordae_tree
    params_Z.chordae         = params.chordae; 
    params_Z.chordae.C_left  = rand(size(params.chordae.C_left)); 
    params_Z.chordae.C_right = rand(size(params.chordae.C_right)); 
end 

Z_linearized = linearize_internal_points(Z, params_Z); 


fprintf('eps\t | taylor series remainder\n'); 


for i = 1:length(epsilon_vals)
    
    ep = epsilon_vals(i); 
    
    % make a new structure for the perturbation 
    params_perturbation = pack_params(X + ep*Z,alpha,beta,N,p_0,R,ref_frac,chordae); 
    
    if chordae_tree 
        params_perturbation.chordae.C_left = params.chordae.C_left + ep*params_Z.chordae.C_left; 
        params_perturbation.chordae.C_right = params.chordae.C_right + ep*params_Z.chordae.C_right; 
    end 
    
    % eval the difference eqns on the perturbation 
    F_perturbed = difference_equations(params_perturbation, filter_params); 
    F_perturbed_linearized = linearize_internal_points(F_perturbed, params); 

    errors(i) = norm(F_perturbed_linearized - F_linearized - ep*J*Z_linearized, 2); 
    
    fprintf('%e\t | %e \n', ep, errors(i)); 

end 

fprintf('\n\n\n\n'); 

figure; 
loglog(epsilon_vals, errors, '-*'); 
hold on 
loglog(epsilon_vals, epsilon_vals.^2, '--'); 

legend('error', 'eps^2')


figure; 

% leaflet part 
for k=1:params.N
    for j=1:params.N
            % in the triangle?
            if (j+k) < (params.N+2)
                
                j
                k
                
                errors = zeros(size(epsilon_vals)); 
                
                
                for i = 1:length(epsilon_vals)

                    ep = epsilon_vals(i); 

                    % make a new structure for the perturbation 
                    params_perturbation = pack_params(X + ep*Z,alpha,beta,N,p_0,R,ref_frac,chordae);
                    
                    if chordae_tree 
                        params_perturbation.chordae.C_left  = params.chordae.C_left  + ep*params_Z.chordae.C_left; 
                        params_perturbation.chordae.C_right = params.chordae.C_right + ep*params_Z.chordae.C_right; 
                    end 

                    % eval the difference eqns on the perturbation 
                    F_perturbed = difference_equations(params_perturbation, filter_params); 
                    F_perturbed_linearized = linearize_internal_points(F_perturbed, params); 

                    diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 
                      
                    range = linear_index_offset(j,k,N) + (1:3);                     
                    errors(i) = norm(diffs(range)); 

                    fprintf('%e\t | %e \n', ep, errors(i)); 

                end 

                fprintf('\n\n\n\n'); 

                loglog(epsilon_vals, errors, '-*'); 
                hold on 
                loglog(epsilon_vals, epsilon_vals.^2, '--'); 

                legend('error', 'eps^2')

            end 
    end 
end 


% chordae part if included 
if chordae_tree 
    
    total_internal = 3*N*(N+1)/2; 
    
    for left_side = [true false]
        for i=1:N_chordae

            left_side
            i

            errors = zeros(size(epsilon_vals)); 


            for ep_idx = 1:length(epsilon_vals)

                ep = epsilon_vals(ep_idx); 

                % make a new structure for the perturbation 
                params_perturbation = pack_params(X + ep*Z,alpha,beta,N,p_0,R,ref_frac,chordae);

                if chordae_tree 
                    params_perturbation.chordae.C_left  = params.chordae.C_left  + ep*params_Z.chordae.C_left; 
                    params_perturbation.chordae.C_right = params.chordae.C_right + ep*params_Z.chordae.C_right; 
                end 

                % eval the difference eqns on the perturbation 
                F_perturbed = difference_equations(params_perturbation, filter_params); 
                F_perturbed_linearized = linearize_internal_points(F_perturbed, params); 

                diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 

                range = range_chordae(total_internal, N_chordae, i, left_side); 
                errors(ep_idx) = norm(diffs(range)); 

                fprintf('%e\t | %e \n', ep, errors(ep_idx)); 

            end 

            fprintf('\n\n\n\n'); 

            loglog(epsilon_vals, errors, '-*'); 
            hold on 
            loglog(epsilon_vals, epsilon_vals.^2, '--'); 

            legend('error', 'eps^2')


        end 
    end 
end 













