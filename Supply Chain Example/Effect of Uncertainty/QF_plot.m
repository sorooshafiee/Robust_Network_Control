clc
clearvars
addpath('../plot tools/');
close all
rng(1000)

T = 24;
num_agents = 3;
P = 1;
K = 4;
F = (2 * rand(P, K) - 1) / K;
B = cell(num_agents,1);
B{1} = 0.5;
B{2} = 0.75;
B{3} = 1;
I_ini = cell(num_agents,1); 
I_ini{1} = 0*ones(P,1);
I_ini{2} = 0*ones(P,1);
I_ini{3} = 0*ones(P,1);

param.P = P;
param.K = K;
param.F = F;
param.B = B;
param.c_B = 1;
param.c_H = 1;
param.I_ini = I_ini;
param.theta_p = 0.2;
param.delay = 0;
theta_d = [0.25, 0.5, 1];
b_l = cell(3,1);
b_u = cell(3,1);

for i = 1 : length(theta_d)
    param.theta_d = theta_d(i);
    [~, b_l{i}, b_u{i}] = local_control(T, P, num_agents, param);
end

plt.font_size = 20;
plt.title = 'Manufacturer - Supplier';
plt.path = '../../Figs/MS';
b_l_2 = [b_l{1}{2}; b_l{2}{2}; b_l{3}{2}];
b_u_2 = [b_u{1}{2}; b_u{2}{2}; b_u{3}{2}];
plot_contracts(b_l_2, b_u_2, T, plt)

plt.title = 'Retailor - Manufacturer';
plt.path = '../../Figs/RM';
b_l_3 = [b_l{1}{3}; b_l{2}{3}; b_l{3}{3}];
b_u_3 = [b_u{1}{3}; b_u{2}{3}; b_u{3}{3}];
plot_contracts(b_l_3, b_u_3, T, plt)
