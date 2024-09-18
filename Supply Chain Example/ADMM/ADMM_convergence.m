clc
clearvars

addpath('../Plot Tools/')
rng(1000)

T = 20;
num_agents = 3;
P = 1;
K = 4;
F = (2 * rand(P, K) - 1) / K;
param.F = F;
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
param.B = B;
param.c_B = 1;
param.c_H = 1;
param.I_ini = I_ini;
param.solver = 'gurobi';
param.theta_p = 0.1;
param.theta_d = 0.5;

param.rho = 0.01;          % parameter rho in ADMM
param.eta = 1;            % parameter eta in ADMM
param.iter_max = 20;
param.delay = 0;

[model, alpha, beta] = local_control(T, P, num_agents, param);
[all_alpha, all_beta] = ADMM(T, P, num_agents, param);
objective = zeros(param.iter_max, 1);
for iter = 1 : param.iter_max
    fprintf('iter %d\n', iter);
    objective(iter) = evaluate_objective(T, P, num_agents, all_alpha{iter}, all_beta{iter}, param);
end
%%
fig = figure;
semilogy(abs(objective - model.get), 'linewidth', 3);
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');
xlabel('\# iterations', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\sum_{i=1}^3 J_i(\beta_i^{(k)}) - J_i(\beta_i^{\star})$', 'Interpreter', 'latex', 'FontSize', 20);
xlim([1,15]);
grid on
saveas(gcf, '../../Figs/ADMM', 'svg')