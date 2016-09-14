% check_jacobian_taylor_series

% checks that the taylor series, using the jacobian included, really works 


epsilon_vals = 10.^(-1:-1:-8); 

errors = zeros(size(epsilon_vals)); 


a = 1; 
r = 1.5;
h = 2; 
N = 9; 
extra = pi/4; 
min_angle = -pi/2 - extra; 
max_angle = -pi/2 + extra; 

filter_params.a = a; 
filter_params.r = r; 
filter_params.h = h;
filter_params.N = N;
filter_params.min_angle = min_angle;
filter_params.max_angle = max_angle;

left = true; 

% reference and initial surfaces are the same 
R = build_reference_surface_commissure(filter_params, left); 

X = R; 
alpha     =  1.0; % spring constants in two directions 
beta      =  1.0;
p_0       = -1.0; 
ref_frac  =  0.5; 

params = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

% difference eqns at center do not change 
F = difference_equations_commissure(params, filter_params, left); 
F_linearized = linearize_internal_points_commissure(F, params); 

% jacobian does not change 
J = build_jacobian_commissure(params, filter_params, left); 

figure; 
spy(J); 
title('Jacobian spy plot for commissural leaflet')


% perturbation also does not change 
Z = zeros(size(X)); 
for j=1:N+2
    for k=1:((N+3)/2)
        if is_internal_commissure(j,k,N)
            Z(:,j,k) = rand(3,1);  
        end 
    end 
end 

params_Z = pack_params(Z,alpha,beta,N,p_0,R,ref_frac); 
Z_linearized = linearize_internal_points_commissure(Z, params_Z); 


fprintf('eps\t | taylor series remainder\n'); 


for i = 1:length(epsilon_vals)
    
    ep = epsilon_vals(i); 
    
    % make a new structure for the perturbation 
    params_perturbation = pack_params(X + ep*Z,alpha,beta,N,p_0,R,ref_frac); 
    
    % eval the difference eqns on the perturbation 
    F_preturbed = difference_equations_commissure(params_perturbation, filter_params, left); 
    F_preturbed_linearized = linearize_internal_points_commissure(F_preturbed, params); 

    errors(i) = norm(F_preturbed_linearized - F_linearized - ep*J*Z_linearized, 2); 
    
    fprintf('%e\t | %e \n', ep, errors(i)); 

end 

fprintf('\n\n\n\n'); 

figure; 
loglog(epsilon_vals, errors, '-*'); 
hold on 
loglog(epsilon_vals, epsilon_vals.^2, '--'); 

legend('error', 'eps^2')


figure; 

for j=1:N+2
    for k=1:((N+3)/2)
        if is_internal_commissure(j,k,N)
                
            j
            k

            errors = zeros(size(epsilon_vals)); 


            for i = 1:length(epsilon_vals)

                ep = epsilon_vals(i); 

                % make a new structure for the perturbation 
                params_perturbation = pack_params(X + ep*Z,alpha,beta,N,p_0,R,ref_frac); 

                % eval the difference eqns on the perturbation 
                F_preturbed = difference_equations_commissure(params_perturbation, filter_params, left); 
                F_preturbed_linearized = linearize_internal_points_commissure(F_preturbed, params); 

                diffs = F_preturbed_linearized - F_linearized - ep*J*Z_linearized; 

                range = linear_index_offset_commissure(j,k,N) + (1:3);                     
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














