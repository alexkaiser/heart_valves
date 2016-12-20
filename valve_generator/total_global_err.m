function err = total_global_err(leaflet)
% 
% Total global error in 2 norm
% 

[F F_chordae_left F_chordae_right] = leaflet.diff_eqns(leaflet); 

if isfield(leaflet, 'repulsive_potential') && leaflet.repulsive_potential
    [F_repulsive F_chordae_left_repulsive F_chordae_right_repulsive] = difference_equations_repulsive(leaflet); 
    F               = F + F_repulsive; 
    F_chordae_left  = F_chordae_left + F_chordae_left_repulsive; 
    F_chordae_right = F_chordae_right + F_chordae_right_repulsive; 
end 


err = norm([F(:); F_chordae_left(:); F_chordae_right(:)]); 

