% check_jacobian_taylor_series


% checks that the taylor series, using the jacobian included, really works 





epsilon_vals = 10.^(-1:-1:-8); 

errors = zeros(size(epsilon_vals)); 


a = 1; 
r = 1.5;
h = 2; 
N = 4; 

filter_params.a = a; 
filter_params.r = r; 
filter_params.h = h;
filter_params.N = N;

% reference and initial surfaces are the same 
R = build_reference_surface(filter_params); 

X = R; 
alpha     =  1.0; % spring constants in two directions 
beta      =  1.0;
p_0       =  1.0; 
ref_frac  =  0.5; 

params = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

% difference eqns at center do not change 
F = difference_equations(params, filter_params); 
F_linearized = linearize_internal_points(F, params); 

% jacobian does not change 
J = build_jacobian(params, filter_params); 

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


params_Z = pack_params(Z,alpha,beta,N,p_0,R,ref_frac); 
Z_linearized = linearize_internal_points(Z, params_Z); 


linear_idx = 1; 
for j=1:params.N
    for k=1:params.N
     
            % in the triangle?
            if (j+k) < (params.N+2)
                
                j
                k
                
                J = build_jacobian_component(params, filter_params,j,k); 
                F = difference_equations_component(params, filter_params,j,k); 
                
                errors = zeros(size(epsilon_vals)); 
                
                for i = 1:length(epsilon_vals)

                    ep = epsilon_vals(i); 

                    % make a new structure for the perturbation 
                    perturbed = X; 
                    perturbed(:,j,k) = perturbed(:,j,k) + eps * Z(:,j,k); 
                    params_perturbation = pack_params(perturbed,alpha,beta,N,p_0,R,ref_frac); 

                    % eval the difference eqns on the perturbation 
                    F_preturbed = difference_equations_component(params_perturbation, filter_params,j,k); 

                    diffs = F_preturbed - F - ep*J*Z(:,j,k); 

                    errors(i) = norm(diffs); 

                    fprintf('%e\t | %e \n', ep, errors(i)); 

                end 

                fprintf('\n\n\n\n'); 

                loglog(epsilon_vals, errors, '-*'); 
                hold on 
                loglog(epsilon_vals, epsilon_vals.^2, '--'); 

                legend('error', 'eps^2')

                linear_idx = linear_idx + 1;

            end 
    end 
end 














