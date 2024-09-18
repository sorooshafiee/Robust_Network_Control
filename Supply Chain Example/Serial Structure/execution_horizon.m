clc
clearvars
addpath('../Plot Tools/')
rng(1000)

P = 1;
K = 4;
T_max = 10;
delay = 0;
num_agents = 12;
theta_p = 0.2;
theta_d = 1;
run_count = 10;
param = cell(run_count, 1);
time_l = zeros(T_max-1, run_count);
time_c = zeros(T_max-1, run_count);

parfor r = 1 : run_count
    fprintf('Running Iteration %d \n', r);
    % initialization for parallel computing
    temp_time_l = zeros(T_max-1, 1);
    temp_time_c = zeros(T_max-1, 1);
    for T = 2 : T_max
        %fprintf('Iteration %d, T = %d \n', r, T);
        F = round(1e3*(2 * rand(P, K) - 1) / K) * 1e-3; % for numerical stability
        B = cell(num_agents,1);
        I_ini = cell(num_agents,1);
        for i = 1 : num_agents
            B{i} = round(1e3 * rand(P, P)) * 1e-3; % for numerical stability
            I_ini{i} = 0 * ones(P,1);
        end
        param{r}.K = K;
        param{r}.delay = delay;
        param{r}.F = F;
        param{r}.B = B;
        param{r}.c_B = round(1e3*rand())*1e-3;  % for numerical stability
        param{r}.c_H = round(1e3*rand())*1e-3;  % for numerical stability
        param{r}.I_ini = I_ini;
        param{r}.delay = delay;
        param{r}.theta_p = theta_p;
        param{r}.theta_d = theta_d;
        % Local model
        t1 = cputime;
        model_l = local_control(T, P, num_agents, param{r});
        t2 = cputime;
        temp_time_l(T-1) =  t2 - t1;
        % Centralized model
        t1 = cputime;
        model_c = centralized_control(T, P, num_agents, param{r});
        t2 = cputime;
        temp_time_c(T-1) = t2 - t1;
    end
    time_l(:,r) = temp_time_l;
    time_c(:,r) = temp_time_c;
end
save execution_horizon T_max time_c time_l
%%
load execution_horizon
plt.prc = 0;
plt.alpha = 0.1;
plt.font_size = 22;
plt.colors = [0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125; 0.50, 0.20, 0.55];
plt.position = [0.35, 0.25, 0.4, 0.3];
plt.grid = true;


plt.lgd = {'Local information', 'Centralized'};
plt.xlabel = 'T (Horizon)';
plt.ylabel = 'Execution time (s)';
plt.path = '../../Figs/execution_horizon';
plot_with_shade(2:T_max, cat(3, time_l, time_c), plt);
