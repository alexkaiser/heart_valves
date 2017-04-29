function T = tension_with_reference(X, X_nbr, R, k_spr, leaflet)
% 
% Returns the tension in the linear constitutive law 
% 
%


if isfield(leaflet, 'collagen_constitutive') && leaflet.collagen_constitutive
    
    collagen_curve       = leaflet.collagen_curve; 
    a                    = collagen_curve.a; 
    b                    = collagen_curve.b; 
    full_recruitment     = collagen_curve.full_recruitment; 
    eta_collagen         = collagen_curve.eta_collagen; 
    collagen_y_intercept = collagen_curve.collagen_y_intercept;
    
    
    E = norm(X - X_nbr)/R - 1.0; 
    
    if E < 0
        T = 0; 
    elseif E < full_recruitment
        T = k_spr * a * (exp(b*E) - 1);
    else 
        T = k_spr * (eta_collagen*E + collagen_y_intercept); 
    end 
    
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