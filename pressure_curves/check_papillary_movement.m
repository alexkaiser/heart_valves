
% 0 

t_diastole_full = .1;  
t_systole_start = .46; % .44;  
t_systole_full  = .5;  

t_cycle_length  = .8; 

times = 0:.001:t_cycle_length; 
vals = zeros(size(times)); 
derivs = zeros(size(times)); 

derivs_linear = zeros(size(times)); 


for i=1:length(times)

    current_time = times(i); 
    
    t_reduced = current_time - t_cycle_length * floor(current_time/(t_cycle_length)); 
    
    
    if (t_reduced < t_diastole_full)
        vals(i)          = (1/2) * (1 - cos(pi * t_reduced/t_diastole_full)); 
        derivs(i)        = pi/(2*t_diastole_full) * sin(pi * t_reduced/t_diastole_full); 
        derivs_linear(i) = 1/t_diastole_full; 
        
        
    elseif  (t_reduced < t_systole_start)
        
        vals(i)          = 1; 
        derivs(i)        = 0; 
        derivs_linear(i) = 0; 
   
    elseif  (t_reduced < t_systole_full)
        
        vals(i)          =  (1/2) * (cos(pi * (t_reduced - t_systole_start)/(t_systole_full - t_systole_start)) + 1); 
        derivs(i)        = -(1/2) *  sin(pi * (t_reduced - t_systole_start)/(t_systole_full - t_systole_start)) * (pi/(t_systole_full - t_systole_start)) ; 
        
        derivs_linear(i) = -1 / (t_systole_full - t_systole_start); 
        
    else
    
        vals(i)          = 0; 
        derivs(i)        = 0;  
        derivs_linear(i) = 0; 
        
    end 

end 


figure; 
plot(times, vals,'k' )
hold on 
plot(times, derivs, 'k--')

plot(times, derivs_linear, 'k-.')

legend('values', 'derivatives', 'piecewise linear deriv')



