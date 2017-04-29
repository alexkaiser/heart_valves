function T_grad = tension_gradient_with_reference(X, X_nbr, R, k_spr, leaflet)
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
    
    E = norm(X - X_nbr)/R - 1.0; 
    
    if E < 0
        coeff = 0; 
    elseif E < full_recruitment
        coeff = k_spr * a * b * exp(b*E);
    else 
        coeff = k_spr * eta_collagen; 
    end 
    
    T_grad = -(coeff/R) * (X_nbr - X) / norm(X_nbr - X);
    
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
    
    