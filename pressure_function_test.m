MMHG_TO_CGS             =  1333.22368
SYSTOLIC_VENTRICULAR_P  =  100.0 % * MMHG_TO_CGS
DIASTOLIC_ATRIAL_P      = -10.0  % * MMHG_TO_CGS

DIASTOLE_TIME = 0.03 
BEAT_TIME = 0.05 

t = 0:.001:.1;
p = zeros(size(t)); 

for i=1:length(t)
    if( (t(i) - BEAT_TIME*floor(t(i)/BEAT_TIME)) < DIASTOLE_TIME )  
        p(i) = DIASTOLIC_ATRIAL_P; 
    else 
        p(i) = SYSTOLIC_VENTRICULAR_P;  
    end 
end 
    

plot(t,p)

