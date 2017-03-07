
p_0 = valve.p_physical / 100; 

valve_linear.anterior.p_0 = p_0; 
valve_linear.posterior.p_0 = p_0; 

if valve_linear.posterior.reflect_pressure
    valve_linear.posterior.p_0 = -p_0; 
end 

[valve_linear.anterior pass_anterior err_anterior] = solve_valve_auto_continuation(valve_linear.anterior, valve_linear.tol_global, valve_linear.max_it, 'anterior'); 

if pass_anterior 
    fprintf('Global solve passed anterior, err = %f\n\n', err_anterior); 
else 
    fprintf('Global solve failed anterior, err = %f\n\n', err_anterior); 
end 

[valve_linear.posterior pass_posterior err_posterior] = solve_valve_auto_continuation(valve_linear.posterior, valve_linear.tol_global, valve_linear.max_it, 'posterior'); 

if pass_posterior 
    fprintf('Global solve passed posterior, err = %f\n\n', err_posterior); 
else 
    fprintf('Global solve failed posterior, err = %f\n\n', err_posterior); 
end