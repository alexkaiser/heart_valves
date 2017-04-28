function T_grad = tension_gradient_with_reference(X, X_nbr, R, k_spr, leaflet)
% 
% Returns the tension in the linear constitutive law 
% 
%


if isfield(leaflet, 'collagen_constitutive') && leaflet.collagen_constitutive
    error('not implemented'); 
    
else 

    % default linear law 
    if ~isfield(leaflet, 'ref_frac')
        ref_frac = leaflet.ref_frac; 
    else 
        ref_frac = 1; 
    end 

    R      = ref_frac * R; 
    
    % Linear law has constant derivative 
    T_grad = -(k_spr/R) * (X_nbr - X) / norm(X_nbr - X);

end 
    
    