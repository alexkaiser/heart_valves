function j_nbr = get_j_nbr(j_nbr_unreduced, k, periodic_j, j_max)
% 
% Returns periodic reduction of j if flag is set at this k 
% Otherwise returns j_nbr_unreduced

j_nbr = j_nbr_unreduced; 

% may have periodic connection in j 
if periodic_j(k)
    if j_nbr == (j_max + 1)
        j_nbr = 1; 
    elseif j_nbr == 0
        j_nbr = j_max; 
    end    
end