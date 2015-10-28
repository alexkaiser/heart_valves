function err = total_global_err(params, filter_params)
% 
% Total global error in 2 norm
% 

F = difference_equations(params, filter_params); 
err = norm(F(:)); 

