clc
clearvars
addpath './plot tools/'
rng(1000)

run_count = 5;
max_agents = 7;
pv_max = 20;
pr = @(t) 18 - tanh(t) + tanh(t - 4) - 2 * tanh (t - 6) + 4 * tanh(t - 17) - 4 * tanh(t - 24);
R = @(t) pv_max/2 * tanh(t - 6) - pv_max/2 * tanh(t - 19);
D_ind = @(t) 100 * sign(t);
D_cmp = @(t) 5 + 10 * tanh(t - 7) - 10 * tanh(t - 18);
D_res = @(t) 5 + 5 * tanh(t - 7) - 5 * tanh(t - 10) + 5 * tanh(t - 17) - 5 * tanh(t - 24);
resolution = 2;
T = 24 / resolution;
B = 60;
sl = -0.5;
penalty = 0.2;

c = pr(1:24); c = sum(reshape(c, [resolution, T]), 1);
pv = R(1:24); pv = sum(reshape(pv, [resolution, T]), 1);
dm_ind = D_ind(1:24); dm_ind = sum(reshape(dm_ind, [resolution, T]), 1);
dm_cmp = D_cmp(1:24); dm_cmp = sum(reshape(dm_cmp, [resolution, T]), 1);
dm_res = D_res(1:24); dm_res = sum(reshape(dm_res, [resolution, T]), 1);
obj_c_s = zeros(max_agents-2, run_count);
obj_c_s_r = zeros(max_agents-2, run_count);
obj_l_s = zeros(max_agents-2, run_count);
obj_l_s_r = zeros(max_agents-2, run_count);
obj_c_c = zeros(max_agents-2, run_count);
obj_c_c_r = zeros(max_agents-2, run_count);
obj_l_c = zeros(max_agents-2, run_count);
obj_l_c_r = zeros(max_agents-2, run_count);
obj_i = zeros(max_agents-2, run_count);
obj_i_r = zeros(max_agents-2, run_count);
c_p = c;

for r = 1 : run_count
    fprintf('Running Iteration %d \n', r);
    temp_obj_c_s = zeros(max_agents-2, 1);
    temp_obj_c_s_r = zeros(max_agents-2, 1);
    temp_obj_l_s = zeros(max_agents-2, 1);
    temp_obj_l_s_r = zeros(max_agents-2, 1);
    temp_obj_c_c = zeros(max_agents-2, 1);
    temp_obj_c_c_r = zeros(max_agents-2, 1);
    temp_obj_l_c = zeros(max_agents-2, 1);
    temp_obj_l_c_r = zeros(max_agents-2, 1);
    temp_obj_i = zeros(max_agents-2, 1);
    temp_obj_i_r = zeros(max_agents-2, 1);
    for num_agents = 3 : max_agents
        fprintf('Iteration %d, # of agents = %d \n', r, num_agents);
        param.pv = pv;
        param.l_pv = -0.2;
        param.u_pv = 0;
        param.dm = [dm_ind; dm_cmp; repmat(dm_res, num_agents - 2, 1)];
        param.l_dm = -0.1;
        param.u_dm = 0.1;
        param.B = B;
        param.c_p = c_p;
        param.c_n = sl * c_p;
        param.I_ini = 0.5 * B * ones(num_agents, 1);
        param.penalty = penalty;
        xi_pv = (param.u_pv - param.l_pv) * rand(num_agents, T) + param.l_pv;
        xi_dm = (param.u_dm - param.l_dm) * rand(num_agents, T) + param.l_dm;
        idx = num_agents - 2;
        
        % Serial
        param.network = diag(ones(num_agents,1), 1) + diag(ones(num_agents,1), -1);
        [temp_obj_c_s(idx), temp_obj_c_s_r(idx)] = rolling_horizon(T, num_agents, param, xi_pv, xi_dm, 'central');
        [temp_obj_l_s(idx), temp_obj_l_s_r(idx)] = rolling_horizon(T, num_agents, param, xi_pv, xi_dm, 'local');
        
        % Fully-connected
        param.network = ones(num_agents, num_agents) - eye(num_agents);
        [temp_obj_c_c(idx), temp_obj_c_c_r(idx)] = rolling_horizon(T, num_agents, param, xi_pv, xi_dm, 'central');
        [temp_obj_l_c(idx), temp_obj_l_c_r(idx)] = rolling_horizon(T, num_agents, param, xi_pv, xi_dm, 'local');
        
        % Individual
        [temp_obj_i(idx), temp_obj_i_r(idx)] = rolling_horizon(T, num_agents, param, xi_pv, xi_dm, 'individual');
    end
    obj_c_s(:,r) = temp_obj_c_s;
    obj_c_s_r(:,r) = temp_obj_c_s_r;
    obj_l_s(:,r) = temp_obj_l_s;
    obj_l_s_r(:,r) = temp_obj_l_s_r;
    obj_c_c(:,r) = temp_obj_c_c;
    obj_c_c_r(:,r) = temp_obj_c_c_r;
    obj_l_c(:,r) = temp_obj_l_c;
    obj_l_c_r(:,r) = temp_obj_l_c_r;
    obj_i(:,r) = temp_obj_i;
    obj_i_r(:,r) = temp_obj_i_r;
end
save results_rolling
%%
load results_rolling

plt.prc = 0;
plt.alpha = 0.1;
plt.font_size = 22;
plt.colors = [0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125; 0.50, 0.20, 0.55];
plt.position = [0.35, 0.25, 0.4, 0.3];
plt.grid = true;
plt.lgd = {'Central (worst-case cost)', 'Local(worst-case cost)', 'Central + Rolling Horizon', 'Local + Rolling Horizon'};
plt.xlabel = '\# of Residential Consumers';
plt.ylabel = 'Cost (\$)';

plt.path = '../Figs/cost_energy_rolling_serial';
plot_with_shade(1:max_agents-2, cat(3, obj_c_s/100, obj_l_s/100, flip(obj_c_s_r)/100, obj_l_s_r/100), plt);
plt.path = '../Figs/cost_energy_rolling_complete';
plot_with_shade(1:max_agents-2, cat(3, obj_c_c/100, obj_l_c/100, flip(obj_c_c_r)/100, flip(obj_l_c_r)/100), plt);