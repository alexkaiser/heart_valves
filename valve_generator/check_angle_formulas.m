

a_ant = 5*pi/6;


a_post = 2*pi - a_ant; 


du = 1e-2; 

u_ant  = 0:du:.5; 

u_post = .5:du:1; 

theta_ant = 2 * a_ant * u_ant - (a_ant/2); 

figure; 
plot(cos(theta_ant), sin(theta_ant)); 
axis equal 
title('anterior, should be right and less than half circle')

theta_post = 2 * a_post * (u_post - .5) - (a_post/2) + pi; 

figure; 

plot(cos(theta_post), sin(theta_post)); 
axis equal 
title('posterior, should be left and little more than half')





