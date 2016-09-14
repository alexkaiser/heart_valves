function err = total_global_err(leaflet)
% 
% Total global error in 2 norm
% 

[F F_chordae_left F_chordae_right] = difference_equations(leaflet); 
err = norm([F(:); F_chordae_left(:); F_chordae_right(:)]); 

