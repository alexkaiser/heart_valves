function T = tension_decreasing(X, X_nbr, ds, c)
% 
% Adjustment to tension from constant 
% 
% Limits to zero at zero and one at infinity 

T = - 1 / (1 +  (c * ds)^(-2) * norm(X_nbr-X)^2); 

