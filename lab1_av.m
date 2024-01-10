% Average velocity experiment
% 

% read dt from URL
url = 'https://raw.githubusercontent.com/tphlabs/data/main/Lab1_AV/test.txt'
block = urlread(url);
C = textscan(block,'%f%f','HeaderLines',1)
S = cell2struct(C,{'n','t'},2)
dt = S.t
%
histfit(dt)
%
N = 100
dt_obs = mean(dt)
dt_obs_err = std(dt) / sqrt(N)
dt_obs_delta = dt_obs_err / dt_obs


%
l = 124.3 *1e-3 % mm
l_err = 0.1*1e-3 % mm
l_delta = l_err / l
s = (402 - 145) * 1e-3 % mm
s_err = 2e-3 % mm
s_delta = s_err / s
L = 127 *1e-2 % cm
L_err = 2 *1e-2 % cm
L_delta = L_err / L
h = 13.8 *1e-3 %mm
h_err = 0.1 *1e-3 %mm
h_delta = h_err / h
g = 9.7949 % m/s^2
% small sin theta = theta
theta = h / L
theta_delta = sqrt(h_delta^2 + L_delta^2)

t2 = sqrt( 2*(s+l) / ( g * theta))
t1 = sqrt( 2*s / ( g * theta))
s_plus_l_delta = (s_err + l_err)/(s+l)
t2_delta = sqrt(s_plus_l_delta^2 + theta_delta^2)/2
t1_delta = (s_delta + theta_delta)/2
t2_err = t2 * t2_delta
t1_err = t1 * t1_delta

dt_model = t2 - t1
dt_model_err = sqrt(t1_err^2 + t2_err^2)
dt_model_delta = dt_model_err / dt_model

AV_model = l / dt_model
AV_model_delta = l_delta + dt_model_delta
AV_model_err = AV_model * AV_model_delta

AV_obs = l / dt_obs
AV_obs_delta = l_delta + dt_obs_delta
AV_obs_err = AV_obs_delta * AV_obs
AV_obs, AV_obs_err

