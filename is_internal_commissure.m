function val = is_internal_commissure(j,k,N)
%  
% Checks whether a given index is an internal point 
% Uses fiber topology of commissural leaflet 
% which is different from anterior or posterior
% 
    if (j < 1) || (k < 1)
        val = false; 
        return 
    elseif j >= k
        val = false; 
        return; 
    elseif ((j+k) >= (N+3))
        val = false;
        return;         
    end 
    
    val = true; 
end 

