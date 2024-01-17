% Instantenious velocity experiment
% IV
clc
clear  
close all
% read counts from URL
url = 'https://raw.githubusercontent.com/tphlabs/data/main/Lab1_IV/d160.txt'

block = urlread(url);
C = textscan(block,'%f%f%f','HeaderLines',1)
S = cell2struct(C,{'n','t','counts'},2)
time = S.t
counts = S.counts
%
plot(counts,'.')
% Uncertainties
% x_er - absolute error (having dimension)
% x_erl - relative error (dimensionless)
% Cart
l = 124.5e-3 % mm
l_er = 0.2e-3 % mm
l_erl = l_er / l
N = 236 
N_er = 1
% Air track
L = 127.5e-2 % cm
L_err = .5e-2 % c,
% Height
h = 16.0e-3 % mm
h_err = 0.2e-3 % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model a = g h / L
g = 9.7949 % m/s^2
g_er = 0.01e-2  % m/s^2
g_err = g_er / g
a0 = - g * h / L
a0_err = sqrt(g_err^2 + h_err^2 + L_err^2)
a0_er = a0 * a0_err
a0_low = a0 - 2 * a0_er
a0_hi = a0 + 2 * a0_er

% segment #1
i0 = 115
i1 = 374
x = counts(i0:i1) * l / N
t = time(i0:i1)


% parabolic regression
figure()
plot(t, x, '.')
[a1, a1_low, a1_high] = parabolicfit(t, x)

% v = dx/dt first derivative
n = 10
[t1, v] = derivative1(t, x, n)

% linear regression
figure()
plot(t1, v)
[fitresult, gof]  = fit(t1, v, 'poly1')
ci = confint(fitresult, 0.67) % 67% confidence = 1 sigma
a2 = fitresult.p1
a2_low = ci(1,1)
a2_hi = ci(2,1)

% a = dx^2/dt^2 second derivative
n = 10
[t1, a] = derivative2(t, x, n)
% linear regression
figure()
plot(t1, a)
a3 = mean(a)
a3_er = std(a)
a3_low = a3 - a3_er
a3_hi = a3 + a3_er

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% up and down comparison
[ xmax, imax] = max(x)
xup = x(imax:-1:1)
tup = t(1:imax) - t(1)
xdown = x(imax:end)
tdown = t(imax:end) - t(imax)
figure()
plot(tup, xup, 'r.')
hold on
plot(tdown, xdown, 'b.')
hold off
% confidence
[aup, aup_low, aup_high] = parabolicfit(tup, xup)
[adown, adown_low, adown_high] = parabolicfit(tdown, xdown)

 



function [t1, v] = derivative1(t, x, n)
    dx = x(1+2*n:end) - x(1:end-2*n)
    dt = t(1+2*n:end) - t(1:end-2*n)
    t1 = t(1+n:end-n)
    v = dx ./ dt
end

function [t1, a] = derivative2(t, x, n)
    dx = x(1+2*n:end) - 2* x(1+n:end-n) + x(1:end-2*n)
    deltat = t(2) - t(1)
    t1 = t(1+n:end-n)
    a = dx ./ (deltat * n)^2
end

function [a, a_low, a_hi] = parabolicfit(t, x)
    [fitresult, gof]  = fit(t, x, 'poly2')
    ci = confint(fitresult, 0.67) % 67% confidence = 1 sigma
    a = fitresult.p1 * 2
    a_low = ci(1,1) * 2
    a_hi = ci(2,1) * 2
end
