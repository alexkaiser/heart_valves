function [] = solve_valve_script(restart_number)


if restart_number ~= 0
    data_str = sprintf('data_iteration_%d', restart_number); 
    load(data_str); 
    start_it = restart_number; 
    max_it_global = 100000; 
    plot_and_save_freq = 1; 
else 

    a = 1; 
    r = 1.5;
    h = 2; 
    N = 32; 

    filter_params.a = a; 
    filter_params.r = r; 
    filter_params.h = h;
    filter_params.N = N;

    % reference and initial surfaces are the same 
    R = build_reference_surface(filter_params); 

    X = R; 
    alpha =  1.0; % spring constants in two directions 
    beta  =  1.0;
    p_0   = -1.0; 

    params = pack_params(X,alpha,beta,N,p_0,R); 

    fig = surf_plot(params, filter_params); 
    title('Reference configuration of surface'); 
    
    'initial difference equation norm'
    err = total_global_err(params, filter_params)
    
    J = build_jacobian(params, filter_params); 
    fig = figure; 
    spy(J); 
    title('jacobian non zero pattern on exact solution')
    
    tol_global = 1e-1; 
    max_it_global = 100000; 
    
    plot_and_save_freq = 1; 
    start_it = 0; 
    
    err_over_time = zeros(max_it_global,1); 
    
end 



random_preturbed_start = false; 
if random_preturbed_start
%     for j=1:N
%         for k=1:N
% 
%             % in the triangle?
%             if (j+k) < (N+2)
%                 X(:,j,k) = X(:,j,k) + 0.001*randn(); 
%             end 
%             
%         end 
%     end 

    X(:,1,1) = X(:,1,1) + 0.1 * randn(); 

    params = pack_params(X,alpha,beta,N,p_0,R); 
end 


J = build_jacobian(params, filter_params); 
fig = figure; 
spy(J); 
title('jacobian non zero pattern on preturbed solution')





[params pass err_over_time it] = solve_valve(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time); 

if pass 
    disp('Global solve passed')
else 
    disp('Global solve failed')
end 


'difference equation norm'
err = total_global_err(params, filter_params)


fig = surf_plot(params, filter_params); 
title('Final time difference equation solution'); 


fig = figure; 
err_over_time = err_over_time(1:it); 
plot(err_over_time); 
title('error through iterations')


