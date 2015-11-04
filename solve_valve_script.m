function [] = solve_valve_script(restart_number)


if restart_number ~= 0
    data_str = sprintf('data_iteration_%d', restart_number); 
    load(data_str); 
    start_it = restart_number; 
    
else 

    a = 1; 
    r = 1.5;
    h = 2; 
    N = 64; 

    filter_params.a = a; 
    filter_params.r = r; 
    filter_params.h = h;
    filter_params.N = N;

    % reference and initial surfaces are the same 
    R = build_reference_surface(filter_params); 

    X = R; 
    alpha     =  1.0; % spring constants in two directions 
    beta      =  1.0;
    p_0       = -20.0; 
    ref_frac  =  0.5; 

    params = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

    fig = surf_plot(params, filter_params); 
    title('Reference configuration of surface'); 
    
    'initial difference equation norm'
    err = total_global_err(params, filter_params)
    
    J = build_jacobian(params, filter_params); 
    fig = figure; 
    spy(J); 
    title('jacobian non zero pattern on initial')
    
    tol_global = 1e-12; 
    max_it_global = 100000; 
    
    plot_and_save_freq = 10; 
    start_it = 0; 
    
    newton_step_coeff = 1.0/256.0; 
    
    err_over_time = zeros(max_it_global,1); 
    
end 


% random_preturbed_start = false; 
% if random_preturbed_start
% %     for j=1:N
% %         for k=1:N
% % 
% %             % in the triangle?
% %             if (j+k) < (N+2)
% %                 X(:,j,k) = X(:,j,k) + 0.001*randn(); 
% %             end 
% %             
% %         end 
% %     end 
% 
%     X(:,1,1) = X(:,1,1) + 0.1 * randn(); 
% 
%     params = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 
% end 
%
% J = build_jacobian(params, filter_params); 
% fig = figure; 
% spy(J); 
% title('jacobian non zero pattern on preturbed solution')



[params pass err_over_time it] = solve_valve(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time, newton_step_coeff); 

if pass 
    disp('Global solve passed')
else 
    disp('Global solve failed')
end 


'difference equation norm'
err = total_global_err(params, filter_params)


fig = surf_plot(params, filter_params); 
title(sprintf('Final time difference equation solution, p = %f, ref_frac = %f' , params.p_0, params.ref_frac )); 
name = sprintf('surf_p_%f', params.p_0); 
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

title(sprintf('error through iterations, p = %f, ref_frac = %f', params.p_0, params.ref_frac))
name = sprintf('error_p_%f', params.p_0)
printfig(fig, strcat(name, '.eps')); 
saveas(fig, strcat(name, '.fig'),'fig'); 

continuation = false; 
if continuation
    
    pressure_vals = [-5.5, -6.0, -6.5, -7.0, -7.5, -8.0]
    
    for p_0 = pressure_vals
    
        params.p_0 = p_0
    
        if pass 

            start_it = 0; 
            
            [params pass err_over_time it] = solve_valve(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time, newton_step_coeff); 

            'difference equation norm'
            err = total_global_err(params, filter_params)


            fig = surf_plot(params, filter_params); 
            
            title(sprintf('Final time difference equation solution, p_0 = %f', p_0)); 
            
            
        else 
            error('cannot run continuation without success ')
        end 
    end 
end 






