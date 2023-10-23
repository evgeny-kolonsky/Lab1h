clc
close all

%  In experimental science, the expression ‘artifact’ is sometimes used to refer 
% to experimental results which are not manifestations of the natural 
% phenomena under investigation, but are due to the particular experimental arrangement.

1. Artifacts, Works, and Authors

N = 240
l = 127 *1e-3 % mm
L = 127 *1e-2 % cm

x = countB / N * l;
t = times1 - times1(1);

fit2 = fit(t, x, 'poly2')
acc2 = fit2.p1 * 2
ci = confint(fit2) % confidence intervals
delta_p1 = (ci(2,1) - ci(1,1))/2
delta_acc2 = delta_p1 *2



% plot with uncertainties
%plot(fit2, 'predobs')
%hold on 
%plot(t,x, '.')

%plot(t, x, '.')

n = 25
dx = x(2*n+1:end) - x(1:end-2*n);
dt = t(2*n+1:end) - t(1:end-2*n);
v = dx./ dt;
t1 = t(n+1:end-n);

%plot(t1,  v, 'r.');


a = g * h / L


%
xl = x(1:249); tl = t(1:249); xr = x(249:end); tr = t(249:end);
fit2l = fit(tl, xl, 'poly2')
fit2r = fit(tr, xr, 'poly2')

plot(t,x, 'b.')
hold on
%plot(fit2l, 'predobs')
plot(fit2r, 'predobs', 0.67)
plot(fit2l, 'predobs', 0.67)
hold off


acc2l = fit2l.p1 * 2
acc2r = fit2r.p1 * 2

%friction
mu = abs(acc2l - acc2r)/2/g




