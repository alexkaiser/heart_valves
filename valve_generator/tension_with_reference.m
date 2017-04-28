function T = tension_with_reference(X, X_nbr, R, k_spr, leaflet)
% 
% Returns the tension in the linear constitutive law 
% 
%


if isfield(leaflet, 'collagen_constitutive') && leaflet.collagen_constitutive
    error('not implemented')
    
else 

    % default linear law 

    if ~isfield(leaflet, 'ref_frac')
        ref_frac = leaflet.ref_frac; 
    else 
        ref_frac = 1; 
    end 

    R      = ref_frac * R; 
    X_norm = norm(X - X_nbr); 

    T = k_spr * (X_norm/R - 1.0); 

end 