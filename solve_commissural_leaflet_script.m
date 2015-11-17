function [] = solve_commissural_leaflet_script(restart_number)


if restart_number ~= 0
    data_str = sprintf('data_iteration_%d', restart_number); 
    load(data_str); 
    start_it = restart_number; 
    
else 

    a = 1; 
    r = 1.5;
    h = 2; 
    N = 4 + 1; 
    extra = pi/2; 
    center = -pi/2; 
    min_angle =  center - extra; 
    max_angle =  center + extra; 
    left = true; 


    filter_params.a = a; 
    filter_params.r = r; 
    filter_params.h = h;
    filter_params.N = N;
    filter_params.min_angle = min_angle;
    filter_params.max_angle = max_angle;


    % reference and initial surfaces are the same 
    R = build_reference_surface_commissure(filter_params, left);

    X = R; 
    alpha =  1.0; % spring constants in two directions 
    beta  =  1.0;
    p_0   = -1.0; % check signs on this 
    ref_frac = 0.8; 

    params_left = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

    fig = surf_plot_commissure(params_left, filter_params, left); 
    title('Reference configuration of surface'); 
    
    'initial difference equation norm'
    err = total_global_err_commissure(params_left, filter_params, left); 
    
    J = build_jacobian_commissure(params_left, filter_params, left); 
    fig = figure; 
    spy(J); 
    title('jacobian non zero pattern on initial commissural leaflet')
    
    tol_global = 1e-12; 
    max_it_global = 100000; 
    
    plot_and_save_freq = 10; 
    start_it = 0; 
    
    err_over_time = zeros(max_it_global,1); 
end 


[params_left pass err_over_time it] = solve_commissure_leaflet(params_left, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time, left); 

if pass 
    disp('Global solve passed')
else 
    disp('Global solve failed')
end 


'difference equation norm'
err = total_global_err_commissure(params_left, filter_params, left)


fig = surf_plot_commissure(params_left, filter_params, left); 
title(sprintf('Final time difference equation solution, p = %f, ref_frac = %f' , params_left.p_0, params_left.ref_frac )); 
name = sprintf('surf_p_%f', params_left.p_0); 
printfig(fig, strcat(name, '.eps'));
saveas(fig, strcat(name, '.fig'),'fig'); 


fig = figure; 
err_over_time = err_over_time(1:it); 
plot(err_over_time); 


fig = figure; 
semilogy(err_over_time); 
hold on 
semilogy(2.^(-(1:it)) , '--'); 

legend('error', '2^-iteration')
title('convergence comparison with 2^{-(iteration)}')
xlabel('iteration')
ylabel('log(err)')

title(sprintf('error through iterations, p = %f, ref_frac = %f', params_left.p_0, params_left.ref_frac))
name = sprintf('error_p_%f', params_left.p_0)
printfig(fig, strcat(name, '.eps')); 
saveas(fig, strcat(name, '.fig'),'fig'); 
