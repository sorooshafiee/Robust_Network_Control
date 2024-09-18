function model = centralized_control(T, P, num_agents, param)

    K = param.K;             % number of factor variables in demand
    F = param.F;             % factors coefficients
    B = param.B;             % production constants
    c_B = param.c_B;         % backlog cost
    c_H = param.c_H;         % holding cost
    I_ini = param.I_ini;     % initial I
    delay = param.delay;     % delay in the supply chain
    theta_p = param.theta_p; % production uncertainty
    theta_d = param.theta_d; % demand uncertainty
    
    % wave function for demand
    wave = zeros(P, T);
    wave(1:2:end, :) = 2 + repmat(cos(2 * pi * (1:T) / (T-1)), ceil(P/2), 1);
    wave(2:2:end, :) = 2 + repmat(sin(2 * pi * (1:T) / (T-1)), fix(P/2), 1);


    % auxiliary function
    vec = @(x) x(:);
    
    % create optimization model
    model = rsome();
    
    % uncertain variables
    xi = cell(num_agents, 1);
    u_set = model.set();
    for i = 1 : num_agents-1
        xi{i} = model.random(P, T);
        u_set.append(xi{i} <= 0);
        u_set.append(xi{i} >= -theta_p);
    end
    xi{num_agents} = model.random(P+K, T);
    xi_d = xi{num_agents}(P+1:end, :);
    u_set.append(xi{num_agents}(1:P,:) <= 0)
    u_set.append(xi{num_agents}(1:P,:) >= -theta_p)
    u_set.append(xi_d <= theta_d)
    u_set.append(xi_d >= -theta_d)
    U = model.ambiguity;
    U.suppset(u_set); 
    model.with(U);
    
    % epigraphic decision variables for objective function
    tau = model.decision(num_agents, 1);
    
    % parameters of the decision rules u and o
    Gamma_xi = cell(num_agents, num_agents, T);
    gamma = cell(num_agents, T);
    Lambda_xi = cell(num_agents, num_agents, T);
    lambda = cell(num_agents, T);
    
    % initialization
    R_t = cell(num_agents, 1);
    I = cell(num_agents, 1);
    u_t = cell(num_agents, 1);
    o_t = cell(num_agents, 1);
    sum_o = cell(num_agents, 1);
    for i = 1 : num_agents
        I{i} = I_ini{i};
        sum_o{i} = 0;
    end

    for t = 1 : T
        
        for i = 1 : num_agents
            
            % decison rule u
            u_t{i} = 0;
            for j = 1 : num_agents
                if t == 1
                    xi_t_j = 0;
                    Gamma_xi{i,j,t} = 0;
                else
                    if j == i
                        xi_t_j = vec(xi{j}(:, 1:t-1));
                        Gamma_xi{i,j,t} = model.decision(P, size(xi_t_j, 1));
                    else
                        if t <= delay + 1
                            xi_t_j = 0;
                            Gamma_xi{i,j,t} = 0;
                        else
                            xi_t_j = vec(xi{j}(:, 1:t-1-delay));
                            Gamma_xi{i,j,t} = model.decision(P, size(xi_t_j, 1));
                        end
                    end
                    
                end
                u_t{i} = u_t{i} + Gamma_xi{i,j,t} * xi_t_j;
            end
            gamma{i,t} = model.decision(P, 1);
            u_t{i} = u_t{i} + gamma{i,t};
            
            % decision rule o
            o_t{i} = 0;
            for j = 1 : num_agents
                xi_t_j = vec(xi{j}(:, 1:t));
                Lambda_xi{i,j,t} = model.decision(P, size(xi_t_j, 1));
                o_t{i} = o_t{i} + Lambda_xi{i,j,t} * xi_t_j;
            end
            lambda{i,t} = model.decision(P, 1);
            o_t{i} = o_t{i} + lambda{i,t};
        end
        
        % add dynamic and other constraints
        for i = 1 : num_agents
            R_t{i} = B{i} * u_t{i} + xi{i}(1:P, t);
            if i < num_agents
                I{i} = I{i} + R_t{i} - u_t{i+1};
            else
                D_t = F * xi_d(:, t) + wave(:, t);
                I{i} = I{i} + R_t{i} - D_t;
            end
            model.append(o_t{i} >= c_H * I{i});
            model.append(o_t{i} >= -c_B * I{i});
            for p = 1 : P
                sum_o{i} = sum_o{i} + o_t{i}(p);
            end
        end
    end
    
    objective = 0;
    for i = 1 : num_agents
        model.append(tau(i) >= sum_o{i});
        objective = objective + tau(i);
    end

    model.min(objective);
    model.Param.solver = 'gurobi';
    model.Param.display = 0;
    model.solve()
end