function objective = evaluate_objective(T, P, num_agents, alpha, beta, param)

    P = param.P;            % number of products
    K = param.K;            % number of factor variables in demand
    F = param.F;            % factors coefficients
    B = param.B;            % production constants
    c_B = param.c_B;        % backlog cost
    c_H = param.c_H;        % holding cost
    I_ini = param.I_ini;    % initial I
    delay = param.delay;    % delay in the supply chain
    
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
    s = cell(num_agents, 1);
    u_set = model.set();
    for i = 1 : num_agents-1
        xi{i} = model.random(P, T);
        s{i+1} = model.random(P, T);
        u_set.append(xi{i} <= 0);
        u_set.append(xi{i} >= -0.1);
        u_set.append(s{i+1} <= 1);
        u_set.append(s{i+1} >= -1);
    end
    xi{num_agents} = model.random(P+K, T);
    xi_d = xi{num_agents}(P+1:end, :);
    u_set.append(xi{num_agents}(1:P,:) <= 0)
    u_set.append(xi{num_agents}(1:P,:) >= -0.1)
    u_set.append(xi_d <= 0.5)
    u_set.append(xi_d >= -0.5)
    U = model.ambiguity;
    U.suppset(u_set); 
    model.with(U);
    
    % epigraphic decision variables for objective function
    tau = model.decision(num_agents, 1);
    
    % parameters of the decision rules u and o
    Gamma_xi = cell(num_agents, T);
    Gamma_s = cell(num_agents, T);
    gamma = cell(num_agents, T);
    Psi_xi = cell(num_agents, T);
    Psi_s = cell(num_agents, T);
    psi = cell(num_agents, T);
    
    % initialization
    xi_t = cell(num_agents, 1);
    xi_t_ = cell(num_agents, 1);
    s_t = cell(num_agents, 1);
    s_t_ = cell(num_agents, 1);
    zeta_t = cell(num_agents, 1);
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
            if t == 1
                xi_t{i} = 0;
                Gamma_xi{i,t} = 0;
            else
                xi_t{i} = vec(xi{i}(:, 1:t-1));
                Gamma_xi{i,t} = model.decision(P, size(xi_t{i}, 1));
            end
            gamma{i,t} = model.decision(P, 1);
            if i < num_agents
                if t <= delay
                    s_t{i+1} = 0;
                    Gamma_s{i+1,t} = 0;
                else
                    s_t{i+1} = vec(s{i+1}(:, 1:t-delay));
                    Gamma_s{i+1,t} = model.decision(P, size(s_t{i+1}, 1));
                end
                u_t{i} = Gamma_xi{i,t} * xi_t{i} + Gamma_s{i+1,t} * s_t{i+1} + gamma{i,t};
            else    
                u_t{i} = Gamma_xi{i,t} * xi_t{i} + gamma{i,t};
            end
            
            % decision rule o
            xi_t_{i} = vec(xi{i}(:, 1:t));
            Psi_xi{i,t} = model.decision(P, size(xi_t_{i}, 1));
            psi{i,t} = model.decision(P, 1);
            if i < num_agents
                s_t_{i+1} = vec(s{i+1}(:, 1:t));
                Psi_s{i+1,t} = model.decision(P, size(s_t_{i+1}, 1));
                o_t{i} = Psi_xi{i,t} * xi_t_{i} + Psi_s{i+1,t} * s_t_{i+1} + psi{i,t};
            else
                o_t{i} = Psi_xi{i,t} * xi_t_{i} + psi{i,t};
            end
        end
        
        % add dynamic and other constraints
        for i = 1 : num_agents
            R_t{i} = B{i} * u_t{i} + xi{i}(1:P, t);
            if i < num_agents
                model.append(u_t{i+1} - beta{i+1}(:, t) <= alpha{i+1}(:,t));
                model.append(beta{i+1}(:, t) - u_t{i+1} <= alpha{i+1}(:,t));
                zeta_t{i+1} = alpha{i+1}(:, t) .* s{i+1}(:, t) + beta{i+1}(:, t);
                I{i} = I{i} + R_t{i} - zeta_t{i+1};
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
    model.Param.solver = param.solver;
    model.Param.display = 0;
    model.solve()
    objective = model.get;
end