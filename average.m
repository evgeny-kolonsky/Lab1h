% Average velocity
% Evgeny Kolonsky Feb 2024

close all
clear

% import measurements to vector t
url = 'https://raw.githubusercontent.com/evgeny-kolonsky/Lab1_Cart/main/cart.txt';
data = readmatrix(url);
t = data(:,2);
figure(1)
hold on
plot(t,'.')

mu = mean(t);
sigma = std(t);
N = length(t);
sigma_mu = sigma / sqrt(N);

% suspected_outlier 
[suspected, ix] = max(abs(t - mu)/sigma/sqrt(2));

P = (1 - cdf('normal', suspected)) + cdf('normal', -suspected)
% 'Chauvene criterion value to be compared with 1/2: {N*P:.1e}')
if N * P < .5
  % 'outlier: to be deleted'
  plot(ix, t(ix), 'rx')
  t(ix) = [];
  % new average and sigma
  mu = mean(t);
  sigma = std(t);
  N = N - 1;
  sigma_mu = sigma / sqrt(N);
end
hold off
legend('data', 'outlier')


figure(2)
txt1 = sprintf('Measured time: \n %.1f +- %.1f ms', mu*1e3, sigma_mu*1e3);
histfit(t * 1e3)
title(txt1)

%% Model
L = 1270e-3;  dL = 1e-3;  %mm
l = 127.7e-3; dl = 1e-3; % mm
s0 = 270e-3; ds0 = 1e-3; %cm
h = 9.1e-3; dh = .1e-3; %mm
g = 9.8; dg = 1e-4; % m/c2

s1 = s0 + l; ds1 = sqrt(dl^2 + ds0^2);


% Relative errors
eL = dL / L;
el = dl / l;
es0 = ds0 / s0;
es1 = ds1 / s1;
eh = dh / h;
eg = dg / g;

% acceleration
a = g * h / L;
ea = sqrt(eg^2 + eh^2 + eL^2);

% time 
t0 = sqrt(2 * s0 /a)
et0 = sqrt(es0^2+ ea^2)/2;
dt0 = t0 * et0;

t1 = sqrt(2 * s1 /a)
et1 = sqrt(es1^2+ ea^2)/2;
dt1 = t1 * et1;

deltat = t1 - t0;
ddeltat = sqrt(dt1^2 + dt0^2);
edeltat = ddeltat /  deltat;

txt2= sprintf('Expected time: \n %.0f +- %.0f s', deltat*1e3, ddeltat*1e3);

