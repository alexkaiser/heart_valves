function p_new = filter_sound(p_0, dt, tau)
% 
% solves ODE 
%    p + tau dp/dt = tau p_0
% 

n_steps = length(p_0) - 1; 
p_new = zeros(length(p_0),1); 
p_new(1) = p_0(1); 

for i=1:n_steps
    p_new(i+1) = (1 - dt/tau) * p_new(i) + (p_0(i+1) - p_0(i)); 
end 

