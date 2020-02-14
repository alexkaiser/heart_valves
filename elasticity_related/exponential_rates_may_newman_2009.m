
% data from "A Hyperelastic Constitutive Law for Aortic Valve Tissue" 
% May-Newman et al, 2009 J Biomech Eng

% Fig 4, bottom left two panels 

% extracted with Engauge Digitizer.app

% stretch, pascals
radial = [1,1.81899e-12
1.462845,163.74
1.504282,264.42
1.527072,385.05
1.544061,493.6
1.552348,577.99
1.562293,742.72
1.575552,983.77
1.581354,1264.94
1.595028,1614.42
1.604972,2084.37
1.614503,2393.67
1.614917,2614.56]; 

% stretch, pascals
circumfernetial = [1.0001558,0.01
1.1341121,108.03
1.1471963,109.08
1.1604361,126.35
1.1750779,204.46
1.1809969,233.28
1.185514,302.47
1.1883178,395.82
1.1956386,679.81
1.1957944,866.06
1.2042056,1150.14
1.2065421,1943.85
1.2113707,2632.5
1.2163551,3369.75
1.223676,4184.1]; 

strain_radial    = radial(:,1) - 1; 
stress_radial_Pa = radial(:,2); 

% manually place data through origin 
strain_radial(1) = 0; 
stress_radial_Pa(1) = 0; 

myfittype = fittype('a*(exp(b*E) - 1.0)',...
    'dependent',{'sigma'},'independent',{'E'},...
    'coefficients',{'a','b'})

myfit_rad = fit(strain_radial,stress_radial_Pa,myfittype)

plot(myfit_rad, strain_radial, stress_radial_Pa)
hold on 

fprintf('Radial exponential rate = %f', myfit_rad.b)


strain_circ      = circumfernetial(:,1) - 1; 
stress_circ_Pa   = circumfernetial(:,2); 
% manually place data through origin 
strain_circ(1) = 0; 
stress_circ_Pa(1) = 0; 

myfittype = fittype('a*(exp(b*E) - 1.0)',...
    'dependent',{'sigma'},'independent',{'E'},...
    'coefficients',{'a','b'})

myfit_circ = fit(strain_circ, stress_circ_Pa, myfittype)

plot(myfit_circ, strain_circ, stress_circ_Pa) 

fprintf('Circ exponential rate = %f\n', myfit_circ.b)

