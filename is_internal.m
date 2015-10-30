function val = is_internal(j,k,N)
%  
% Checks whether a given index is an internal point 
% 

    if (j < 1) || (k < 1)
        val = false; 
        return 
    elseif (j+k) >= (N+2)   
        val = false; 
        return; 
    end 
    
    val = true; 
end 

