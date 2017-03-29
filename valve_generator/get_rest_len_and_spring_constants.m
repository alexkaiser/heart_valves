function [k R] = get_rest_len_and_spring_constants(X, X_nbr, tension, strain, leaflet)
% 
% Given two points, a tension between them, and a desired specified strain 
% Compute rest length and spring constant that gives the current force 
% 
% Input: 
%     X          Point 
%     X_nbr      Other point 
%     tension    Force between the two (not that this must be an actual force, not a density)
%     strain     Desired strain 
% 
% Output: 
%     k          Spring constant k 
%                multiplies RELATIVE strain to get force 
%                This implies that k has units of force 
%     R          Rest length 
% 

if isfield(leaflet, 'collagen_constitutive') && leaflet.collagen_constitutive
    
    if strain < collagen_curve.full_recruitment
        error('Must set constants with affine part of constitutive law for now. Exponential not implemented'); 
    end 
    
    % affine part of constitutive law 
    k = tension / (collagen_curve.eta_collagen * strain + collagen_curve.collagen_y_intercept);
else 
    % linear law by default 
    k = tension / strain; 
end

% Rest length is determined by strain 
R = norm(X - X_nbr) / (strain + 1);

