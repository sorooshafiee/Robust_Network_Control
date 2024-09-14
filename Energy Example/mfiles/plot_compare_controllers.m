clc
clearvars
addpath './plot tools/'
addpath '../csv/'
rng(1000)

T = 24;
resolution = 2;
T = T / resolution;
run_count = 10;
units = [1, 2, 3, 9, 10, 11];
max_agents = length(units);
pv = zeros(max_agents, T);
r_pv = zeros(max_agents, T);
dm = zeros(max_agents, T);
r_dm = zeros(max_agents, T);
B = zeros(max_agents, 1);

for i = 1 : max_agents
    % PV
    fname = strcat('pv_power_unit_', int2str(units(i)), '.csv');
    raw_data = csvread(fname) / 360;
    pv(i, :) = sum(reshape(raw_data(1, :), [resolution, T]), 1);
    tmp = raw_data(2, :) ./ (raw_data(1, :) + 1e-12);
    r_pv(i, :) = mean(reshape(tmp, [resolution, T]), 1);
    
    % Demand
    fname = strcat('load_power_unit_', int2str(units(i)), '.csv');
    raw_data = csvread(fname) / 360;
    dm(i, :) = sum(reshape(raw_data(1, :), [resolution, T]), 1);
    tmp = raw_data(2, :) ./ (raw_data(1, :) + 1e-12);
    r_dm(i, :) = mean(reshape(tmp, [resolution, T]), 1);
    
    % Battery
    fname = strcat('battery_unit_', int2str(units(i)), '.csv');
    raw_data = csvread(fname);
    B(i) = raw_data(1);
end
pr = @(t) 18 - tanh(t) + tanh(t - 4) - 2 * tanh (t - 6) + 4 * tanh(t - 17) - 4 * tanh(t - 24);
c = pr(1:24);
c = sum(reshape(c, [resolution, T]), 1);
sl = -0.5;
penalty = 0.2;

obj_c_s = zeros(max_agents-1, run_count);
time_c_s = zeros(max_agents-1, run_count);
obj_c_c = zeros(max_agents-1, run_count);
time_c_c = zeros(max_agents-1, run_count);
obj_l_s = zeros(max_agents-1, run_count);
time_l_s = zeros(max_agents-1, run_count);
obj_l_c = zeros(max_agents-1, run_count);
time_l_c = zeros(max_agents-1, run_count);
obj_i = zeros(max_agents-1, run_count);
time_i = zeros(max_agents-1, run_count);

for r = 1 : run_count
    fprintf('Running Iteration %d \n', r);
    temp_obj_c_s = zeros(max_agents-1, 1);
    temp_time_c_s = zeros(max_agents-1, 1);
    temp_obj_c_c = zeros(max_agents-1, 1);
    temp_time_c_c = zeros(max_agents-1, 1);
    temp_obj_l_s = zeros(max_agents-1, 1);
    temp_time_l_s = zeros(max_agents-1, 1);
    temp_obj_l_c = zeros(max_agents-1, 1);
    temp_time_l_c = zeros(max_agents-1, 1);
    temp_obj_i = zeros(max_agents-1, 1);
    temp_time_i = zeros(max_agents-1, 1);
    c_p = c .* (0.9 + 0.2 * rand(1,T));
    for num_agents = 3 : max_agents
        fprintf('Iteration %d, # of agents = %d \n', r, num_agents);
        param.pv = pv(1:num_agents, :);
        param.l_pv = -r_pv(1:num_agents, :)/2;
        param.u_pv = r_pv(1:num_agents, :)/2;
        param.dm = dm(1:num_agents, :);
        param.l_dm = -r_dm(1:num_agents, :)/2;
        param.u_dm = r_dm(1:num_agents, :)/2;
        param.B = B;
        param.c_p = c_p;
        param.c_n = sl * c_p;
        param.I_ini = 0.0 * B;
        param.penalty = penalty;
        idx = num_agents - 1;
        
        % Central + Serial Network
        param.network = diag(ones(num_agents,1), 1) + diag(ones(num_agents,1), -1);
        [temp_obj_c_s(idx), ~, ~, ~, temp_time_c_s(idx)] = central_control(T, num_agents, param);
        
        % Central + Complete Network
        param.network = ones(num_agents, num_agents) - eye(num_agents);
        [temp_obj_c_c(idx), ~, ~, ~, temp_time_c_c(idx)] = central_control(T, num_agents, param);
        
        % Local + Serial Network
        param.network = diag(ones(num_agents,1), 1) + diag(ones(num_agents,1), -1);
        [temp_obj_l_s(idx), ~, ~, ~, temp_time_l_s(idx)] = local_control(T, num_agents, param);
        
        % Local + Complete Network
        param.network = ones(num_agents, num_agents) - eye(num_agents);
        [temp_obj_l_c(idx), ~, ~, ~, temp_time_l_c(idx)] = local_control(T, num_agents, param);
        
        % Individual
        [temp_obj_i(idx), ~, ~, temp_time_i(idx)] = individual_control(T, num_agents, param);
    end
    obj_c_s(:,r) = temp_obj_c_s;
    time_c_s(:,r) = temp_time_c_s;
    obj_c_c(:,r) = temp_obj_c_c;
    time_c_c(:,r) = temp_time_c_c;
    obj_l_s(:,r) = temp_obj_l_s;
    time_l_s(:,r) = temp_time_l_s;
    obj_l_c(:,r) = temp_obj_l_c;
    time_l_c(:,r) = temp_time_l_c;
    obj_i(:,r) = temp_obj_i;
    time_i(:,r) = temp_time_i;
end
save results_controller
%%
% load results_controller
% 
plt.prc = 0;
plt.alpha = 0.0;
plt.font_size = 22;
plt.colors = [0, 0.45, 0.75; 0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125];
plt.position = [0.35, 0.25, 0.4, 0.3];
plt.grid = true;
plt.lgd = {'Central + Serial Network', 'Central + Fully-Connected Network', 'Local + Serial Network', 'Local + Fully-Connected Network', 'Individual'};
plt.xlabel = '\# of Residential Consumers';
plt.ylabel = 'Execution time (s)';
plt.path = '../Figs/execution_energy';
plot_with_shade_modified(1:max_agents-2, cat(3, time_c_s, time_c_c, time_l_s, time_l_c, time_i), plt);
plt.ylabel = 'Optimal Objective Value (\$)';
plt.path = '../Figs/cost_energy';
plot_with_shade_modified(1:max_agents-2, cat(3, obj_c_s, obj_c_c, obj_l_s, obj_l_c, obj_i)/100, plt);

%%
% load results_controller
plt.prc = 0;
plt.alpha = 0.1;
plt.font_size = 22;
plt.colors = [0, 0.45, 0.75; 0.85, 0.325, 0.01; 0.925, 0.70, 0.125];
plt.position = [0.35, 0.25, 0.4, 0.3];
plt.grid = true;
plt.lgd = {'Central', 'Local', 'Individual'};
plt.xlabel = '\# of Residential Consumers';
plt.ylabel = 'Optimization time (s)';

plt.path = '../Figs/execution_energy_serial';
plot_with_shade(1:max_agents-2, cat(3, time_c_s, time_l_s, time_i), plt);

plt.path = '../Figs/execution_energy_complete';
plot_with_shade(1:max_agents-2, cat(3, time_c_c, time_l_c, time_i), plt);

plt.ylabel = 'Worst-case cost (\$)';
plt.path = '../Figs/cost_energy_serial';
plot_with_shade(1:max_agents-2, cat(3, obj_c_s, obj_l_s, obj_i)/100, plt);

plt.path = '../Figs/cost_energy_complete';
plot_with_shade(1:max_agents-2, cat(3, obj_c_c, obj_l_c, obj_i)/100, plt);