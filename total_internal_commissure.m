function tot = total_internal_commissure(N)
% 
% Returns the total number of internal points 
% in the commissural leaflet topology 
% 

if mod(N,2) ~= 1
    error('N must be odd for commisural leaflets'); 
end 

tot = ((N+1)/2)^2 + 1; 

