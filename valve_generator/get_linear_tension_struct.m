function Q = get_linear_tension_struct(X, X_nbr, R, R_nbr, kappa, ref_frac, j, k, j_nbr, k_nbr)
%
% Returns a struct for tension with the following fields
% 
%     val        Tension 
%     j          Current index 
%     k          Current index 
%     j_nbr      Neighboring index 
%     k_nbr      Neighboring index 
%     G          Jacobian, which is a gradient since tensions are scalar valued 
% 
% Input: 
%     X          Jacobian is taken with respect to this variable 
%     X_nbr      Relevant neighbor in X
%     R          Reference coordinate at current location 
%     R_nbr      Reference coordinate at neighbor location 
%     kappa      Spring constant 
%     ref_frac   Multiplier for rest length 
%     j          Current index 
%     k          Current index 
%     j_nbr      Neighboring index 
%     k_nbr      Neighboring index
% 
% 
% 

X_norm =            norm(X - X_nbr);
R_norm = ref_frac * norm(R - R_nbr);

Q.val = kappa * (X_norm/R_norm - 1.0); 

Q.j = j; 
Q.k = k; 
Q.j_nbr = j_nbr; 
Q.k_nbr = k_nbr; 

Q.G = (kappa/R_norm) * (X - X_nbr) / X_norm; 

