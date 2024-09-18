clc
clearvars
addpath('../Plot Tools/')

P = 1;
K = 4;
T = 5;
num_agents_max = 12;
run_count = 10;
theta_p = 0.2;
theta_d = 1;
param = cell(run_count, 1);
obj_l = zeros(num_agents_max-1, run_count);
obj_c = zeros(num_agents_max-1, run_count);
obj_l_d = zeros(num_agents_max-1, run_count);
obj_c_d = zeros(num_agents_max-1, run_count);

parfor r = 1 : run_count
    fprintf('Running Iteration %d \n', r);
    % initialization for parallel computing
    temp_obj_l = zeros(num_agents_max-1, 1);
    temp_obj_c = zeros(num_agents_max-1, 1);
    temp_obj_l_d = zeros(num_agents_max-1, 1);
    temp_obj_c_d = zeros(num_agents_max-1, 1);
    for num_agents = 2 : num_agents_max
        % fprintf('Iteration %d, # of agents = %d \n', r, num_agents);
        F = round(1e3*(2 * rand(P, K) - 1) / K) * 1e-3; % for numerical stability
        B = cell(num_agents,1);
        I_ini = cell(num_agents,1);
        for i = 1 : num_agents
            B{i} = 1 - 0.25 * round(1e3 * rand(P, P)) * 1e-3; % for numerical stability
            I_ini{i} = 0 * ones(P,1);
        end
        param{r}.K = K;
        param{r}.F = F;
        param{r}.B = B;
        param{r}.c_B = round(1e3*rand())*1e-3;  % for numerical stability
        param{r}.c_H = round(1e3*rand())*1e-3;  % for numerical stability
        param{r}.I_ini = I_ini;
        param{r}.solver = 'gurobi';
        param{r}.theta_p = theta_p;
        param{r}.theta_d = theta_d;
        % Local model
        param{r}.delay = 0;
        model_l = local_control(T, P, num_agents, param{r});
        temp_obj_l(num_agents-1) = model_l.get;
        param{r}.delay = 1;
        model_l_d = local_control(T, P, num_agents, param{r});
        temp_obj_l_d(num_agents-1) = model_l_d.get;
        % Centralized model
        param{r}.delay = 0;
        model_c = centralized_control(T, P, num_agents, param{r});
        temp_obj_c(num_agents-1) = model_c.get;
        param{r}.delay = 1;
        model_c_d = centralized_control(T, P, num_agents, param{r});
        temp_obj_c_d(num_agents-1) = model_c_d.get;
    end
    obj_l(:, r) = temp_obj_l;
    obj_l_d(:,r) = temp_obj_l_d;
    obj_c(:, r) = temp_obj_c;
    obj_c_d(:,r) = temp_obj_c_d;
end
save suboptimality_serial_agents num_agents_max obj_l obj_l_d obj_c obj_c_d
%%
load suboptimality_serial_agents
plt.prc = 0;
plt.alpha = 0.1;
plt.font_size = 26;
plt.colors = [0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125; 0.50, 0.20, 0.55];
plt.position = [0.35, 0.25, 0.4, 0.3];
plt.grid = true;

plt.lgd = {'Without Delay', 'With Delay'};
plt.xlabel = 'N (Manufacturer)';
plt.ylabel = 'Suboptimality (\%)';
plt.path = '../../Figs/suboptimality_serial_agents';
plot_with_shade(0:num_agents_max-2, ...
                cat(3, abs((obj_l - obj_c) ./ obj_c) * 100, ...
                abs((obj_l_d - obj_c_d) ./ obj_c_d) * 100), plt);
