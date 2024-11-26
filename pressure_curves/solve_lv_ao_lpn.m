function [times, P_lv, Q_ao, P_ao, V_lv] = solve_lv_ao_lpn(dt, t_final, V0_ventricle, Vrd, Vrs, Emax, Emin, ... 
                                              Q_in, act, P_ao_initial, R_proximal, C, R_distal, R_av)


n_times = t_final/dt; 

% state variable arrays
times = zeros(n_times, 1); 
P_lv  = zeros(n_times, 1); 
Q_ao  = zeros(n_times, 1); 
P_ao  = zeros(n_times, 1); 
V_lv  = zeros(n_times, 1); 

P_wk  = zeros(n_times, 1); 
Vrest = zeros(n_times, 1); 
Elas  = zeros(n_times, 1); 

% initial conditions 
times(1) = 0; 
P_lv(1)  = 0; 
Q_ao(1)  = 0; 
P_ao(1)  = P_ao_initial;
V_lv(1)  = V0_ventricle; 

P_wk(1)  = P_ao_initial;
Vrest(1) = 0; 
Elas(1)  = 0;

for n = 1:(n_times-1)

    times(n+1) = dt*n;
    
    % aortic valve is simple diode 
    if P_lv(n) > P_ao(n)
        Q_ao(n+1) = (1/R_av) * (P_lv(n) - P_ao(n)); 
    else 
        Q_ao(n+1) = 0; 
    end 
    
    % rcr update 
    P_wk(n+1) = P_wk(n) + (dt/C) * (-P_wk(n)/R_distal + Q_ao(n)); 
    P_ao(n+1) = P_wk(n+1) + R_proximal * Q_ao(n); 

    % eqn 3 for ventricular volume 
    V_lv(n+1) = V_lv(n) + dt * (Q_in(times(n)) - Q_ao(n)); 

%     if t < 0.2
%         act_temp = 0; 
%     else 
%         act_temp = act(n+1); 
%     end 

    act_temp = act(times(n));

    Vrest(n+1) = (1.0 - act_temp) * (Vrd - Vrs) + Vrs;
    Elas(n+1) = (Emax - Emin) * act_temp + Emin;

    % update pressure with new volumes 
    P_lv(n+1) = Elas(n+1) * (V_lv(n+1) - Vrest(n+1));     
end 







