function val = in_domain_with_bc(j,k,N)
%  
% Checks whether a given index is a valid point
% Includes boundary values 
% Uses fiber topology of commissural leaflet 
% which is different from anterior or posterior
% 
    if (j < 1) || (k < 1)
        val = false; 
        return 
    elseif k > j
        val = false; 
        return; 
    elseif ((j+k) > (N+3))
        val = false;
        return;         
    end 
    
    val = true; 
end 

