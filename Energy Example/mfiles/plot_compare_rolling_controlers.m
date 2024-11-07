clc
clearvars
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
        idx = num_agents - 2;
        xi_pv = (param.u_pv - param.l_pv) .* rand(num_agents, T) + param.l_pv;
        xi_dm = (param.u_dm - param.l_dm) .* rand(num_agents, T) + param.l_dm;
        
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

figure('units','normalized','outerposition',[0 0 0.8 0.4]) 
labels_1 = {'CS','LS','D'};
labels_2 = {'M = 3','M = 4','M = 5','M = 6'}; 
colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880];
x = {[obj_c_s_r, obj_c_s(:,1)]', [obj_l_s_r, obj_l_s(:,1)]', [obj_i_r, obj_i(:,1)]'};

h = boxplotGroup(x, 'PrimaryLabels',labels_1, ...
'SecondaryLabels',labels_2,'GroupLabelType','Vertical','Widths',0.75, ...
'Colors',colors,'GroupType','betweenGroups'); 
set(gca,'FontSize', 22);
set(gca, 'TickLabelInterpreter', 'latex');
set(findobj(h.boxplotGroup,'tag','Outliers'),'LineWidth',3,'Marker','*','MarkerSize',22)
set(findobj(h.boxplotGroup,'tag','Upper Whisker'),'LineWidth',3)
set(findobj(h.boxplotGroup,'tag','Lower Whisker'),'LineWidth',3)
set(findobj(h.boxplotGroup,'tag','Upper Adjacent Value'),'LineWidth',3)
set(findobj(h.boxplotGroup,'tag','Lower Adjacent Value'),'LineWidth',3)
set(findobj(h.boxplotGroup,'tag','Box'),'LineWidth',3)
set(findobj(h.boxplotGroup,'tag','Median'),'LineWidth',3)
ylabel('Cost (\$)', 'Interpreter', 'latex');
saveas(gcf, './Figs/compare_rolling', 'svg')


