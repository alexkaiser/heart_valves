function [times, P_lv, Q_ao, P_ao, V_lv, R_valve_series] = solve_lv_ao_lpn(dt, t_final, V0_ventricle, Vrd, Vrs, Emax, Emin, ... 
                                              Q_in, act, P_ao_initial, R_proximal, C, R_distal, R_av, R_av_closed, steepness_av)


n_times = floor(t_final/dt); 

% state variable arrays
times = zeros(n_times, 1); 
P_lv  = zeros(n_times, 1); 
Q_ao  = zeros(n_times, 1); 
P_ao  = zeros(n_times, 1); 
V_lv  = zeros(n_times, 1); 

P_wk  = zeros(n_times, 1); 
Vrest = zeros(n_times, 1); 
Elas  = zeros(n_times, 1); 

R_valve_series = zeros(n_times, 1); 

% initial conditions 
times(1) = 0; 
P_lv(1)  = 0; 
Q_ao(1)  = 0; 
P_ao(1)  = P_ao_initial;
V_lv(1)  = V0_ventricle; 

P_wk(1)  = P_ao_initial;
Vrest(1) = 0; 
Elas(1)  = 0;

tanh_valve = true; 

linear_ramp_valve = false; 
valve_resistance_ramp_time = 0.05; 

r_valve_increment = (R_av_closed - R_av) * (dt/valve_resistance_ramp_time);

r_valve_temp = R_av_closed;

% valve_open = false; 
% valve_closed = true; 
% valve_opening = false; 
% valve_closing = false; 

bump_radius = valve_resistance_ramp_time / 2;
h = (2/pi) * bump_radius; 
cos_bump = @(x) (abs(x) <= bump_radius) .* (1/h) .* (2/pi) .* cos(x/h).^2; 
mesh_increments = linspace(-bump_radius,bump_radius,valve_resistance_ramp_time/dt);
increments = cos_bump(mesh_increments);


valve_state = 1; 
% closed, opening, open, closing 

inc_idx = 1;


for n = 1:(n_times-1)

    times(n+1) = dt*n;
 

    
    if tanh_valve
        delta_p_valve = P_lv(n) - P_ao(n); 
        R_tanh_valve = R_av + (R_av_closed - R_av) * 0.5 * (1 + tanh(steepness_av * -delta_p_valve));
        R_valve_series(n) = R_tanh_valve;
        Q_ao(n+1) = (1/R_tanh_valve) * delta_p_valve;
    elseif linear_ramp_valve
        
        delta_p_valve = P_lv(n) - P_ao(n); 

        if valve_state == 1 
            % closed 
            r_valve_temp = R_av_closed; 
        elseif valve_state == 2
            % opening 
            r_valve_temp = r_valve_temp - (R_av_closed - R_av) * dt*increments(inc_idx);
            inc_idx = inc_idx + 1;
        elseif valve_state == 3
            % open
            r_valve_temp = R_av;
        elseif valve_state == 4
            % closing 
            r_valve_temp = r_valve_temp + (R_av_closed - R_av) * dt*increments(inc_idx)
            inc_idx = inc_idx + 1
        end 
        
        % check for changes in state 
        if (valve_state == 1) && (delta_p_valve > 0)
            % if closed and positive pressure difference, start opening 
            times(n+1)
            valve_state = 2
            inc_idx = 1; 
        end 
        
        if (valve_state == 2) && (inc_idx > length(increments))
            % state open, found opening complete 
            times(n+1)
            valve_state = 3
        end 
        
        if (valve_state == 3) && (delta_p_valve < 0) 
            % if open and negative pressure differnce, start closing 
            times(n+1)
            valve_state = 4
            inc_idx = 1;
        end 
        
        if (valve_state == 4) && (inc_idx > length(increments))
            % state closed, found closing complete
            times(n+1)
            valve_state = 1
        end

        % linear ramp valve, generally bad 
%         if delta_p_valve > 0
%             % resistance is decreasing 
%             if (r_valve_temp > R_av) && (inc_idx <= length(increments))
%                 r_valve_temp = r_valve_temp - (R_av_closed - R_av) * dt*increments(inc_idx);
%                 inc_idx = inc_idx + 1;
%             else
%                 % steady state 
%                 r_valve_temp = R_av;
%                 
%                 % increments start anew 
%                 inc_idx = 1; 
%             end 
% 
%         end
%         
%         if delta_p_valve < 0
%             
%             if n>4798
%                 'pause'
%             end 
%             
%             % resistance is decreasing 
%             if (r_valve_temp < R_av_closed) && (inc_idx <= length(increments))
%                 r_valve_temp = r_valve_temp + (R_av_closed - R_av) * dt*increments(inc_idx);
%                 inc_idx = inc_idx + 1
%             else
%                 % steady state 
%                 r_valve_temp = R_av_closed;
%                 
%                 % increments start anew 
%                 inc_idx = 1; 
%             end 
%         end     
        
        R_valve_series(n) = r_valve_temp; 
        
        Q_ao(n+1) = (1/r_valve_temp) * (P_lv(n) - P_ao(n)); 
        
    else 
        % aortic valve is simple diode 
        if P_lv(n) > P_ao(n)
            Q_ao(n+1) = (1/R_av) * (P_lv(n) - P_ao(n)); 
        else 
            Q_ao(n+1) = 0; 
        end 
    end 
    
    % rcr update 
    back_euler_rcr = true; 
    if back_euler_rcr
        P_wk(n+1) = ((C/dt) * P_wk(n) + Q_ao(n+1)) / (C/dt + 1.0/R_distal);        
    else 
        P_wk(n+1) = P_wk(n) + (dt/C) * (-P_wk(n)/R_distal + Q_ao(n)); 
    end
    P_ao(n+1) = P_wk(n+1) + R_proximal * Q_ao(n+1);

    % eqn 3 for ventricular volume 
    V_lv(n+1) = V_lv(n) + dt * (Q_in(times(n+1)) - Q_ao(n+1)); 

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







