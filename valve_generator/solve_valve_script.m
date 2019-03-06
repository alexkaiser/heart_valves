function [] = solve_valve_script(restart_number)

% Copyright (c) 2019, Alexander D. Kaiser
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


if restart_number ~= 0
    data_str = sprintf('data_iteration_%d', restart_number); 
    load(data_str); 
    start_it = restart_number; 
    
else 

    a = 1; 
    r = 1.5;
    h = 2; 
    N = 64; 
    extra = .9*pi/4; 
    min_angle = -pi/2 - extra; 
    max_angle =  pi/2 + extra; 

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
    p_0       = -30.0; 
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
    
    newton_step_coeff = 1.0; %1.0/256.0; 
    
    err_over_time = zeros(max_it_global,1); 
    
end 


[params pass err_over_time it] = solve_valve(params, filter_params, tol_global, max_it_global, plot_and_save_freq, start_it, err_over_time); 

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
semilogy(err_over_time, '*-'); 
hold on 
title('Convergence of Newton iteration')
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






