function Q = get_prescribed_tension_struct(kappa, j, k, j_nbr, k_nbr)
%
% Returns a struct for tension with the following fields
% 
%     val        Tension 
%     j          Current index 
%     k          Current index 
%     j_nbr      Neighboring index 
%     k_nbr      Neighboring index 
% 
% Input: 
%     X          Jacobian is taken with respect to this variable 
%     X_nbr      Relevant neighbor in X
%     kappa      Spring constant 
%     j          Current index 
%     k          Current index 
%     j_nbr      Neighboring index 
%     k_nbr      Neighboring index
% 

Q.val = kappa; 
Q.j = j; 
Q.k = k; 
Q.j_nbr = j_nbr; 
Q.k_nbr = k_nbr; 


