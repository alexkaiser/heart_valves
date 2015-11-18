function [] = solve_valve_script_four_leaflet(restart_number)


if restart_number ~= 0
    data_str = sprintf('data_iteration_%d', restart_number); 
    load(data_str); 
    start_it = restart_number; 
    
else 

    % symmetric no overlap geometry for now 
    commissural_angle = pi/6; 
    
    
    a = 1; 
    r = 1.5;
    h = 2; 
    N = 32; 
    min_angle_posterior = -(pi/2 - commissural_angle/2); 
    max_angle_posterior =   pi/2 - commissural_angle/2; 

    filter_params_posterior.a = a; 
    filter_params_posterior.r = r; 
    filter_params_posterior.h = h;
    filter_params_posterior.N = N;
    filter_params_posterior.min_angle = min_angle_posterior;
    filter_params_posterior.max_angle = max_angle_posterior;

    % reference and initial surfaces are the same 
    R = build_reference_surface(filter_params_posterior); 

    X = R; 
    alpha     =  1.0; % spring constants in two directions 
    beta      =  1.0;
    p_0       = -20.0; 
    ref_frac  =  0.5; 

    params_posterior = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

    fig = surf_plot(params_posterior, filter_params_posterior); 
    title('Reference configuration of posterior surface'); 
    
    'initial difference equation norm'
    err = total_global_err(params_posterior, filter_params_posterior)
    
    tol_global = 1e-12; 
    max_it_global = 100000; 
    
    plot_and_save_freq = 10; 
    start_it = 0; 
    
    newton_step_coeff = 1.0; %1.0/256.0; 
    
    err_over_time_posterior = zeros(max_it_global,1); 
    
    
    % anterior 
    less = .2*pi/4; % shrink by this much 
    min_angle_anterior = -pi/2 + commissural_angle/2; 
    max_angle_anterior =  pi/2 + commissural_angle/2; 

    filter_params_anterior.a = a; 
    filter_params_anterior.r = r; 
    filter_params_anterior.h = h;
    filter_params_anterior.N = N;
    filter_params_anterior.min_angle = min_angle_anterior;
    filter_params_anterior.max_angle = max_angle_anterior;

    % reference and initial surfaces are the same 
    R = build_reference_surface(filter_params_anterior); 
    X = R; 

    params_anterior = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

    fig = surf_plot(params_anterior, filter_params_anterior); 
    title('Reference configuration of anterior surface'); 
    
    'initial difference equation norm'
    err = total_global_err(params_anterior, filter_params_anterior)
    
    err_over_time_anterior = zeros(max_it_global,1); 
end 




[params_posterior pass err_over_time_posterior it] = solve_valve(params_posterior, filter_params_posterior, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time_posterior); 

% reflect posterior to actually have two leaflets
params_posterior.X(1,:,:) = -params_posterior.X(1,:,:); 

if pass 
    disp('Global solve passed posterior')
else 
    disp('Global solve failed')
end


[params_anterior pass err_over_time_anterior it] = solve_valve(params_anterior, filter_params_anterior, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time_anterior, newton_step_coeff); 

if pass 
    disp('Global solve passed anterior')
else 
    disp('Global solve failed')
end 


'difference equation norm'
err = total_global_err(params_posterior, filter_params_posterior)

fig = figure; 
fig = surf_plot(params_posterior, filter_params_posterior, fig);
hold on 
fig = surf_plot(params_anterior, filter_params_anterior, fig);
title(sprintf('Final time difference equation solution, p = %f, ref_frac = %f' , params_posterior.p_0, params_posterior.ref_frac )); 
name = sprintf('surf_p_%f', params_posterior.p_0); 
printfig(fig, strcat(name, '.eps'));
saveas(fig, strcat(name, '.fig'),'fig'); 



fig = figure; 
semilogy(err_over_time_posterior, '*-'); 
hold on 
semilogy(err_over_time_anterior, 's-'); 

legend('error posterior', 'error anterior', 'location', 'SouthWest'); 
xlabel('iteration')
ylabel('log(err)')

title(sprintf('error through iterations, p = %f, ref frac = %f', params_posterior.p_0, params_posterior.ref_frac))
name = sprintf('error_p_%f', params_posterior.p_0)
printfig(fig, strcat(name, '.eps')); 
saveas(fig, strcat(name, '.fig'),'fig'); 



% continuation = false; 
% if continuation
%     
%     pressure_vals = [-5.5, -6.0, -6.5, -7.0, -7.5, -8.0]
%     
%     for p_0 = pressure_vals
%     
%         params.p_0 = p_0
%     
%         if pass 
% 
%             start_it = 0; 
%             
%             [params pass err_over_time it] = solve_valve(params, filter_params_posterior, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time, newton_step_coeff); 
% 
%             'difference equation norm'
%             err = total_global_err(params, filter_params_posterior)
% 
% 
%             fig = surf_plot(params, filter_params_posterior); 
%             
%             title(sprintf('Final time difference equation solution, p_0 = %f', p_0)); 
%             
%             
%         else 
%             error('cannot run continuation without success ')
%         end 
%     end 
% end 
% 





