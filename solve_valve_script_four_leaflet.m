function [] = solve_valve_script_four_leaflet(restart_number)


if restart_number ~= 0
    data_str = sprintf('data_iteration_%d', restart_number); 
    load(data_str); 
    start_it = restart_number; 
    
else 

    % symmetric for now 
    commissure_angle = 0.0; 
    
    % note that if extra posterior is equal to commissure angle
    % then the poster leaflet has angle pi 
    extra_posterior = pi/6; 
    
    a = 1; 
    r = 1.5;
    h = 2; 
    N = 32; 

    min_angle_posterior = -(pi/2 - commissure_angle/2 + extra_posterior/2); 
    max_angle_posterior =  (pi/2 - commissure_angle/2 + extra_posterior/2); 

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
    p_0       = -10.0; 
    ref_frac  =  0.5; 

    params_posterior = pack_params(X_posterior,alpha,beta,N,p_0,R_posterior,ref_frac); 

    fig = surf_plot(params_posterior, filter_params_posterior); 
    title('Reference configuration of posterior surface'); 
    
    'initial difference equation norm'
    err = total_global_err(params_posterior, filter_params_posterior)
    
    tol_global = 1e-11; 
    max_it_global = 10; 
    
    plot_and_save_freq = 10; 
    start_it = 0; 
    
    err_over_time_posterior = zeros(max_it_global,1); 
    
    
    % anterior 
    min_angle_anterior = -(pi/2 - commissure_angle/2 - extra_posterior/2); 
    max_angle_anterior =  (pi/2 - commissure_angle/2 - extra_posterior/2); 

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
    
end 


% p_range = -(0:2.5:50); 
% ref_frac_range = .1:.1:1; 

p_range = -(0:2.5:5); 
ref_frac_range = .5:.1:.6; 


for ref_frac = ref_frac_range
    
    % reset the whole thing when the reference fraction changes
    params_posterior = pack_params(X_posterior,alpha,beta,N,p_0,R_posterior,ref_frac); 
    params_anterior = pack_params(X_anterior,alpha,beta,N,p_0,R_anterior,ref_frac); 
    if commissure_angle > 0.0
        params_left = pack_params(X_left,alpha,beta,N_left,p_0,R_left,ref_frac); 
        params_right = pack_params(X_right,alpha,beta,N_right,p_0,R_right,ref_frac); 
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

        if commissure_angle > 0.0
            left = true; 
            fig = surf_plot_commissure(params_left, filter_params_left, left, fig); 

            left = false;
            fig = surf_plot_commissure(params_right, filter_params_right, left, fig); 
        end 
        
        title(sprintf('Final time difference equation solution, p = %f, ref_frac = %f' , params_posterior.p_0, params_posterior.ref_frac )); 
        
        % side view 
        view([0 0]);
        name = sprintf('surf_side_p_%f_ref_frac_%f', p_0, ref_frac); 
        printfig(fig, strcat(name, '.eps'));
        saveas(fig, strcat(name, '.fig'),'fig'); 

        % top view 
        view([0 90]);
        name = sprintf('surf_top_p_%f_ref_frac_%f', p_0, ref_frac); 
        printfig(fig, strcat(name, '.eps'));
        saveas(fig, strcat(name, '.fig'),'fig');         

        % error 
        fig = figure; 
        semilogy(err_over_time_posterior, '*-'); 
        hold on 
        semilogy(err_over_time_anterior, 's-'); 

        if commissure_angle > 0.0
            semilogy(err_over_time_left, 'o-'); 
            semilogy(err_over_time_right, '+-'); 
        end 

        if commissure_angle > 0.0
            legend('posterior', 'anterior', 'left', 'right', 'location', 'SouthWest'); 
        else
            legend('posterior', 'anterior', 'location', 'SouthWest'); 
        end 

        xlabel('iteration')
        ylabel('log(err)')

        title(sprintf('error through iterations, p = %f, ref frac = %f', params_posterior.p_0, params_posterior.ref_frac))
        name = sprintf('error_p_%f_ref_frac_%f', p_0, ref_frac)
        printfig(fig, strcat(name, '.eps')); 
        saveas(fig, strcat(name, '.fig'),'fig'); 

        
        % re-reflect posterior to be back where is needed to run the thing 
        params_posterior.X(1,:,:) = -params_posterior.X(1,:,:); 
        
        close all 
    end
end






