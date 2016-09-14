function err = total_global_err_commissure(params, filter_params, left)
% 
% Total global error in 2 norm
% 

F = difference_equations_commissure(params, filter_params, left); 
err = norm(F(:)); 

