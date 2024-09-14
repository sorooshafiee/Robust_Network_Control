function [obj_w, obj_a] = rolling_horizon(T, num_agents, param, xi_pv, xi_dm, controller)
    I_t = param.I_ini;
    obj_a = 0;
    for t =  0 : T-1
        if strcmp(controller, 'central')
            [obj, g_p, g_n, u] = central_control(T-t, num_agents, param);
        elseif strcmp(controller, 'local')
            [obj, g_p, g_n, u] = local_control(T-t, num_agents, param);
        else
            [obj, g_p, g_n] = individual_control(T-t, num_agents, param);
            u = zeros(num_agents, num_agents);
        end
        g_p = get_values(g_p,false);
        g_n = get_values(g_n,false);
        u = get_values(u,true);
        if t == 0
            obj_w = obj;
        end    
        for i = 1 : num_agents            
            temp_i = g_p(i) - g_n(i);
            for j = 1 : num_agents
                temp_i = temp_i + u(j,i) - u(i,j);
            end
            I_t(i) = I_t(i) + temp_i + param.pv(i, 1) * (1 + xi_pv(i, t+1)) - param.dm(i, 1) * (1  + xi_dm(i, t+1));
            obj_a = obj_a + param.c_p(1) * g_p(i) - param.c_n(1) * g_n(i);
        end
        param.I_ini = I_t;
        param.pv(:, 1) = [];
        param.l_pv(:, 1) = [];
        param.u_pv(:, 1) = [];
        param.dm(:, 1) = [];
        param.l_dm(:, 1) = [];
        param.u_dm(:, 1) = [];
        param.c_p(1) = [];
        param.c_n(1) = [];
    end
end

function val = get_values(var, is_2d)
   if isnumeric(var)
       val = var;
   else
       num_agents = size(var, 1);
       if ~is_2d
           val = zeros(num_agents, 1);
           for i = 1 : num_agents
               if ~isnumeric(var{i,1})
                   val(i) = var{i,1}.get;
               end
           end
       else
           val = zeros(num_agents, num_agents);
           for i = 1 : num_agents
               for j = 1 : num_agents
                   if ~isnumeric(var{i,j,1})
                       val(i,j) = var{i,j}.get;
                   end
               end
           end
       end
   end
end