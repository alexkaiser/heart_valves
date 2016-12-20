function F = diff_eqns_and_linearize(leaflet, handle)


[F_2d F_chordae_left F_chordae_right] = handle(leaflet); 

F = linearize_internal_points(leaflet, F_2d, F_chordae_left, F_chordae_right);

