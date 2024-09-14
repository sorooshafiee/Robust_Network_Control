close all
clear all
%
flagVar = 1;
Np = 10;
Nv = 5;
x = sdpvar(Nv,Np,'full');
x0 = sdpvar(Nv,1);
u = sdpvar(Nv,Np,'full');
w = sdpvar(Nv,Np,'full');
wE = zeros(Nv,Np);
p1 = sdpvar;
p2 = sdpvar;
% create optimization problem
obj = 0;
constr = [];
constrW = [];
for kk = 1:Np
    obj = obj + 3*ones(1,Nv)*u(:,kk);
    if kk == 1
        constr = constr + [x(:,kk) == 0.9*x0 + u(:,kk) + w(:,kk),...
               u(:,kk) >= -1,[u(:,kk) <= 1]];
    else
        constr = constr + [-5<= x(:,kk), x(:,kk) <= 5,... 
               x(:,kk) == 0.9*x(:,kk-1) + u(:,kk) + w(:,kk),...
               u(:,kk) >= -1,u(:,kk) <= 1];
    end
    constrW = constrW + [p2<= w(:,kk), w(:,kk)<=p1];
end

% worst case formulation
tau = sdpvar(1,1);
constr =  constr + [[obj <= tau]:'dualTest'];
obj = tau;

%%  test case
p1 = RobustProblem(obj,constr,constrW,[x0,x],u,[],[],w,wE,[x0;p1;p2]);
p1.robustify(flagVar,0);
p1.getRobust([1*ones(Nv,1);1.5;-1.5]);
p1.solveRobust('gurobi',0);
p1.scenario;
p1.getScenarioLB;
p1.solveScenarioLB('gurobi',0);


price = p1.getDualConstraint('dualTest');
%% Previous impelementation
% Robustify the problem
[optRobust1,realVarRobust1] = getRobust(constr,obj,u,[],[],x,...
    w,wE,1.5*ones(size(w)),-1.5*ones(size(w)),0,...
    flagVar,x0);

% formulate and solve the optimization problem
% replace the parameters to get the correct ordering
paramVar = replace(optRobust1.paramVar,x0,1);

% create the linear optimization problem
conRobustEq = optRobust1.Aeq*paramVar == optRobust1.beq;
conRobustIn = optRobust1.Ain*paramVar <= optRobust1.bin;

objRobust = replace(optRobust1.objRobust,x0,1);

optimize(conRobustEq+conRobustIn,objRobust,sdpsettings('solver','gurobi','verbose',0));
fprintf('Centralized Problem: Old method: %.6f\n',value(objRobust));

