function J_tension = tension_jacobian(X_current,X_nbr,R_current,R_nbr,k_spr,ref_frac)
% 
% Computes the contribution to the Jacobian for X and its neighbor
% 
% Input: 
%     X_current      Jacobian is taken with respect to this variable 
%     X_nbr          Relevant neighbor in X
%     R_current      Reference coordinate at current location 
%     R_nbr          Reference coordinate at neighbor location 
%     k              Spring constant 
% 
% Output: 
%     J_tension      3x3 Jacobian for this tension term 
%                    Signs are NOT included, local will get a negative 
% 

%     J_tension = zeros(3,3);
% 
%     % u difference norms which are used repeatedly
%     X_norm = norm(X_nbr - X_current);
%     R_norm = ref_frac * norm(R_nbr - R_current);
% 
%    
%     for l=1:3
%         for m=1:3
% 
%             J_tension(l,m) = (X_nbr(l) - X_current(l)) * (X_nbr(m) - X_current(m)) / (X_norm^3) ; 
%                     
%             % diagonal term has an extra term
%             if l == m
%                 J_tension(l,m) = J_tension(l,m) + (1.0/R_norm - 1.0/X_norm);                       
%             end
%         end
%     end
%     
%     J_tension = k_spr * J_tension;

    X_norm = norm(X_nbr - X_current);
%    R_norm = ref_frac * norm(R_nbr - R_current);
    
    % jacobian is an outer product 
    % plus a multiple of the identity 
    J_tension = -k_spr / (X_norm^3) * (X_nbr - X_current)*((X_nbr - X_current)') ... 
                -tension_linear(X_current,X_nbr,R_current,R_nbr,k_spr,ref_frac) * eye(3); 
     
end 

