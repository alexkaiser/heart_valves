function [] = solve_valve_script(restart_number)


if restart_number ~= 0
    data_str = sprintf('data_iteration_%d', restart_number); 
    load(data_str); 
    start_it = restart_number; 
    max_it_global = 100000; 
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
    p_0   = -1.0; % check signs on this 

    params = pack_params(X,alpha,beta,N,p_0,R); 

    fig = surf_plot(params, filter_params); 
    title('Reference configuration of surface'); 
    
    'initial difference equation norm'
    err = total_global_err(params, triangle_domain, linear_constitutive)
    
    tol_global = 1e-1; 
    max_it_global = 100000; 
    
    
    
    
end 



[params pass err_over_time] = solve_valve(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time)


if pass 
    disp('Global solve passed')
else 
    disp('Global solve failed')
end 


'difference equation norm'
err = total_global_err() 


fig = surf_plot(params, filter_params); 
title('Final time difference equation solution'); 


fig = figure; 
plot(err_over_time); 
title('error through iterations')


