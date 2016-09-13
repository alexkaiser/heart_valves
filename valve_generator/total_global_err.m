function err = total_global_err(params, filter_params)
% 
% Total global error in 2 norm
% 

[F F_chordae_left F_chordae_right] = difference_equations(params, filter_params); 
err = norm([F(:); F_chordae_left(:); F_chordae_right(:)]); 

