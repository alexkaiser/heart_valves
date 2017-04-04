function err = total_global_err(leaflet)
% 
% Total global error in 2 norm
% 

F = leaflet.diff_eqns(leaflet); 

err = norm(F); 

