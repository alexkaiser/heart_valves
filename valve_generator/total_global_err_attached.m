function err = total_global_err_attached(valve)
% 
% Total global error in 2 norm
% 

[F_anterior F_posterior F_chordae_left F_chordae_right] = difference_equations_bead_slip_attached(valve); 
err = norm([F_anterior(:); F_posterior(:); F_chordae_left(:); F_chordae_right(:)]); 

