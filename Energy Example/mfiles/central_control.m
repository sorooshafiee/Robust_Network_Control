function [objective, g_p, g_n, u, solution_time] = central_control(T, num_agents, param)
    network = param.network;
    pv = param.pv;
    l_pv = param.l_pv;
    u_pv = param.u_pv;
    dm = param.dm;
    l_dm = param.l_dm;
    u_dm = param.u_dm;
    B = param.B;
    c_p = param.c_p;
    c_n = param.c_n;
    I_ini = param.I_ini;
    penalty = param.penalty;
    
    % create optimization model
    model = rsome();
    
    % uncertain variables
    xi_set = model.set();
    xi_pv = cell(num_agents, 1);
    xi_dm = cell(num_agents, 1);
    for i = 1 : num_agents
        xi_pv{i} = model.random(1, T);
        xi_set.append(xi_pv{i} >= l_pv(i,:));
        xi_set.append(xi_pv{i} <= u_pv(i,:));
        xi_dm{i} = model.random(1, T);
        xi_set.append(xi_dm{i} >= l_dm(i,:));
        xi_set.append(xi_dm{i} <= u_dm(i,:));
    end
    XI = model.ambiguity;
    XI.suppset(xi_set); 
    model.with(XI);
    
    % initialization
    I_t = cell(num_agents, 1);
    u = cell(num_agents, num_agents, T);
    g_p = cell(num_agents, T);
    g_n = cell(num_agents, T);
    sum_g = cell(num_agents, 1);
    for i = 1 : num_agents
        I_t{i} = I_ini(i);
        sum_g{i} = 0;
    end
    
    % Let's begin now
    for t = 1 : T
        for i = 1 : num_agents
            % decison rules g_p and g_n
            g_p{i,t} = model.decision(1);
            g_n{i,t} = model.decision(1);
            if t > 1
                for j = 1 : num_agents
                    g_p{i,t}.affadapt(xi_pv{j}(:, 1:t-1));
                    g_n{i,t}.affadapt(xi_pv{j}(:, 1:t-1));
                    g_p{i,t}.affadapt(xi_dm{j}(:, 1:t-1));
                    g_n{i,t}.affadapt(xi_dm{j}(:, 1:t-1));
                end
            end
            % decision rule u
            for j = 1 : num_agents
                if network(i,j) == 1
                    u{j,i,t} = model.decision(1);
                    if t > 1
                        for k = 1 : num_agents
                            u{j,i,t}.affadapt(xi_pv{k}(:, 1:t-1));
                            u{j,i,t}.affadapt(xi_dm{k}(:, 1:t-1));
                        end
                    end
                else
                    u{j,i,t} = 0;
                end
            end
        end
        
        % add constraints on u
        for i = 1 : num_agents
            for j = 1 : num_agents
                if network(i,j) == 1
                    model.append(u{j,i,t} >= 0);
                end
            end
        end
        
        % add dynamic and other constraints
        for i = 1 : num_agents
            temp_i = g_p{i,t} - g_n{i,t};
            for j = 1 : num_agents
                if network(i,j) == 1
                    temp_i = temp_i + u{j,i,t} - u{i,j,t};
                end
            end
            I_t{i} = I_t{i} + temp_i + pv(i, t) * (1 + xi_pv{i}(:, t)) - dm(i, t) * (1  + xi_dm{i}(:, t));
            model.append(I_t{i} >= 0);
            model.append(I_t{i} <= B(i));
            model.append(g_p{i,t} >= 0);
            model.append(g_n{i,t} >= 0);
            sum_g{i} = sum_g{i} + c_p(t) * g_p{i,t} - c_n(t) * g_n{i,t};
            for j = 1 : num_agents
                sum_g{i} = sum_g{i} + penalty * c_p(t) * u{j,i,t};
            end
        end
    end
    
    % epigraphic decision variables for objective function
    tau = model.decision(num_agents, 1);
    objective = 0;
    for i = 1 : num_agents
        model.append(tau(i) >= sum_g{i});
        objective = objective + tau(i);
    end

    model.min(objective);
    model.Param.solver = 'gurobi';
    model.Param.display = 0;
    model.solve()
    objective = model.get;
    solution_time = model.Solution.time;
    clearvars -except objective g_p g_n u solution_time
end