function [leaflet pass err E] = solve_valve_optimization(leaflet, tol)
% 
% Solves valve using built in matlab optimization routine 
% Leaflet must include an energy, difference equations and Jacobian 
% Note that difference equations are the gradient of energy
% and that the Jacobian is the Hessian of energy 
% 
% 


energy_gradient_hessian_handle = @(X) energy_gradient_hessian(X, leaflet); 

% current configuration is initial guess 
X_linearized = linearize_internal_points(leaflet, leaflet.X, leaflet.chordae.C_left, leaflet.chordae.C_right); 

% should probably set some options here... 
options = optimset('LargeScale',        'on', ...
                   'Diagnostics',       'on', ...
                   'Display',           'iter-detailed', ...    % 'iter-detailed'
                   'GradObj',           'on', ...
                   'TolFun',             tol, ...
                   'TolX',               1.0e-12, ...
                   'Hessian',           'on'); 


[X_soln, E, exitflag, output, grad, hessian] = fminunc(energy_gradient_hessian_handle, X_linearized, options); 

'algorithm output:'
output
'algorithm output message:'
output.message


err = norm(grad); 

fprintf('flag = %d,\t E = %f,\t || F || = %f\n', exitflag, E, err); 


% send data back 
leaflet = internal_points_to_2d(X_soln, leaflet); 


% E_copy = leaflet.energy(leaflet); 
% err = total_global_err(leaflet); 
% 
% fprintf('E_copy = %f,\t || F || = %f\n', E_copy, err); 

pass = (err < tol); 









