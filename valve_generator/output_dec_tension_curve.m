function [] = output_dec_tension_curve()

dx = 1e-2; 
x = 0:dx:5; 

f = @(x) 1 - (1 + (x.^2)).^(-1); 

fig = figure;
plot(x, f(x), 'k')

set(gcf,'color',[1 1 1])

xlabel('x')
ylabel('f(x)')

printfig(fig, 'dec_tension_curve.eps')
