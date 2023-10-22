clc
close all



N = 240
l = 127 *1e-3 % mm
L = 127 *1e-2 % cm

x = countB / N * l;
t = times1;

%plot(t, x, '.')

n = 25
dx = x(2*n+1:end) - x(1:end-2*n);
dt = t(2*n+1:end) - t(1:end-2*n);
v = dx./ dt;
t1 = t(n+1:end-n);

plot(t1,  v, 'r.');


a = g * h / L


