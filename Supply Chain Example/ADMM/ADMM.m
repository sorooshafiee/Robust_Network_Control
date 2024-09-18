function [all_alpha, all_beta] = ADMM(T, P, num_agents, param)

    rho = param.rho;            % parameter rho in ADMM
    eta = param.eta;            % parameter eta in ADMM
    iter_max = param.iter_max;  % maximum number of iterations in ADMM
    
    % Initialize the algorithm
    c = 2 * ones(num_agents); c(1) = 1; c(num_agents) = 1;
    lambda = cell(num_agents, 1);
    alpha = cell(num_agents, 1);
    beta = cell(num_agents, 1);
    for i = 1 : num_agents
        lambda{i} = zeros(2*c(i)*P,T);
        alpha{i} = zeros(P,T);
        beta{i} = zeros(P,T);
    end
    z_tilde = exchange_information(alpha, beta);
    
    all_alpha = cell(iter_max,1);
    all_beta = cell(iter_max,1);
    for iter = 1 : iter_max
        fprintf('Iteration %d\n', iter)
        param.lambda = lambda;     
        param.z_tilde = z_tilde;
        % step 1: solve subproblems
        x = cell(num_agents, 1);
        for i = 1 : num_agents
            x{i} = ADMM_subproblem(T, P, num_agents, i, param);
        end
        % step 2: update proximal points
        [alpha, beta] = combine_information(x);
        z_tilde = exchange_information(alpha, beta);
        % step 3: update dual variables
        for i = 1 : num_agents
            lambda{i} = lambda{i} + eta * rho * (x{i} - z_tilde{i});
        end
        % save the results
        all_alpha{iter} = alpha;
        all_beta{iter} = beta;
    end
    
end

function x_i = ADMM_subproblem(T, P, num_agents, i, param)

    K = param.K;                % number of factor variables in demand
    F = param.F;                % factors coefficients
    B = param.B;                % production constants
    c_B = param.c_B;            % backlog cost
    c_H = param.c_H;            % holding cost
    I_ini = param.I_ini;        % initial I
    delay = param.delay;        % delay in the supply chain
    lambda = param.lambda;      % dual variables
    z_tilde = param.z_tilde;    % proximal points
    rho = param.rho;            % parameter rho in ADMM
    
    % wave function for demand
    wave = zeros(P, T);
    wave(1:2:end, :) = 2 + repmat(cos(2 * pi * (1:T) / (T-1)), ceil(P/2), 1);
    wave(2:2:end, :) = 2 + repmat(sin(2 * pi * (1:T) / (T-1)), fix(P/2), 1);

    % auxiliary function
    vec = @(x) x(:);
    
    % create optimization model
    model = rsome();
    
    % uncertain variables
    u_set = model.set();
    if i < num_agents
        xi = model.random(P, T);
        s = model.random(P, T);
        u_set.append(xi <= 0);
        u_set.append(xi >= -0.1);
        u_set.append(s <= 1);
        u_set.append(s >= -1);        
    else
        xi = model.random(P+K, T);
        u_set.append(xi(1:P,:) <= 0)
        u_set.append(xi(1:P,:) >= -0.1)
        u_set.append(xi(P+1:end,:) <= 0.5)
        u_set.append(xi(P+1:end,:) >= -0.5)
    end
    U = model.ambiguity;
    U.suppset(u_set); 
    model.with(U);
    
    % epigraphic decision variables for objective function
    tau = model.decision(1);
    r = model.decision(1);
    
    % decision variables in subproblems
    c = 2 * ones(num_agents); c(1) = 1; c(num_agents) = 1;
    x_i = model.decision(2*c(i)*P, T);
    
    % initialization
    I = I_ini{i};
    sum_o = 0;

    for t = 1 : T
        
        % decison rule u
        gamma = model.decision(P, 1);
        if t == 1
            xi_t = 0;
            Gamma_xi = 0;
        else
            xi_t = vec(xi(:, 1:t-1));
            Gamma_xi = model.decision(P, size(xi_t, 1));
        end
        if i < num_agents
            if t <= delay
                s_t = 0;
                Gamma_s = 0;
            else
                s_t = vec(s(:, 1:t-delay));
                Gamma_s = model.decision(P, size(s_t, 1));
            end
            u_t = gamma + Gamma_xi * xi_t + Gamma_s * s_t;
        else    
            u_t = gamma + Gamma_xi * xi_t;
        end
            
        % decision rule o
        xi_t_ = vec(xi(:, 1:t));
        psi = model.decision(P, 1);
        Psi_xi = model.decision(P, size(xi_t_, 1));
        if i < num_agents
            s_t_ = vec(s(:, 1:t));
            Psi_s = model.decision(P, size(s_t_, 1));
            o_t = psi + Psi_xi * xi_t_ + Psi_s * s_t_;
        else
            o_t = psi + Psi_xi * xi_t_;
        end
        
        % add dynamic and other constraints
        R_t = B{i} * u_t + xi(1:P, t);
        if i == 1
            % simple renaming for clarity
            alpha_i_next = x_i(1:P, :);
            beta_i_next = x_i(P+1:end, :);
            
            zeta_t = alpha_i_next(:, t) .* s(:, t) + beta_i_next(:, t);
            I = I + R_t - zeta_t;
        elseif i == num_agents
            % simple renaming for clarity
            alpha_i = x_i(1:P, :);
            beta_i = x_i(P+1:end, :);
            
            D_t = F * xi(P+1:end, t) + wave(:, t);
            I = I + R_t - D_t;
            model.append(u_t - beta_i(:, t) <= alpha_i(:,t));
            model.append(beta_i(:, t) - u_t <= alpha_i(:,t));
        else
            % simple renaming for clarity
            alpha_i = x_i(1:P, :);
            beta_i = x_i(P+1:2*P, :);
            alpha_i_next= x_i(2*P+1:3*P, :);
            beta_i_next = x_i(3*P+1:end, :);
            
            zeta_t = alpha_i_next(:, t) .* s(:, t) + beta_i_next(:, t);
            I = I + R_t - zeta_t;
            model.append(u_t - beta_i(:, t) <= alpha_i(:,t));
            model.append(beta_i(:, t) - u_t <= alpha_i(:,t));
        end
        
        model.append(o_t >= c_H * I);
        model.append(o_t >= -c_B * I);
        for p = 1 : P
            sum_o = sum_o + o_t(p);
        end
    end
    model.append(tau >= sum_o);
    model.append(r >= vec(lambda{i})' * vec(x_i) + ...
                      0.5 * rho * sumsqr(vec(x_i - z_tilde{i})));
    objective = tau + r;
    model.min(objective);
    model.Param.solver = param.solver;
    model.Param.display = 0;
    model.solve()
    
    x_i = x_i.get;
end

function [alpha, beta] = combine_information(x)
    num_agents = size(x, 1);
    P = size(x{1}, 1) / 2;
    T = size(x{1}, 2);
    x{1} = [nan(2*P, T); x{1}];
    x{num_agents} = [x{num_agents}; nan(2*P, T)];
    alpha = cell(num_agents, 1);
    beta = cell(num_agents, 1);
    for i = 2 :num_agents
        alpha{i} = 0.5 * x{i-1}(2*P+1:3*P,:) + 0.5 * x{i}(1:P,:);
        beta{i} = 0.5 * x{i-1}(3*P+1:4*P,:) + 0.5 * x{i}(P+1:2*P,:);
    end
end

function z_tilde = exchange_information(alpha, beta)
    num_agents = size(alpha, 1);
	z_tilde = cell(num_agents, 1);
    z_tilde{1} = [alpha{2}; beta{2}];
    z_tilde{num_agents} = [alpha{num_agents}; beta{num_agents}];
    for i = 2 : num_agents-1
        z_tilde{i} = [alpha{i}; alpha{i+1}; ...
                      beta{i}; beta{i+1}];
    end
end