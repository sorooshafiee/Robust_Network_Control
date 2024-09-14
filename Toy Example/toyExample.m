clc
close all
clearvars

addpath(genpath('./Robust_Optimization/'))
addpath('../Plot Tools/');

fDual = 0;
fPlot = 1;

xu = 8;
xl = -8;
uu = 4;
ul = -4;

% System matrices
B = [1;0.8];
C = [1,-1;-1,1];
D = [1,0;0,-2];

O1 = [1,-1];
O2 = [1,-1];


% Agent 1
x1 = sdpvar(2,2,'full');
u1 = sdpvar(1,1);
w1 = sdpvar(2,1);

conE1 = [x1(:,2) == B*u1 + C*w1];
conX1 = [xl <= x1(:,2) <= xu];
conU1 = [ul <= u1, u1 <= uu];

con1 = conE1 + conX1 + conU1;
conW1 = [-1<=w1, w1<= 1];

obj1 = O1*x1(:,2);

tau1 = sdpvar(1,1);
objP1 = tau1;
con1 = con1 + [obj1 <= tau1];


% Agent 2
x2 = sdpvar(2,2,'full');
u2 = sdpvar(1,1);
w2 = sdpvar(2,1);

conE2 = [x2(:,2) == B*u2 + C*w2 + D*x1(:,2)];
conX2 = [xl <= x2(:,2) <= xu];
conU2 = [ul <= u2, u2 <= uu];

con2 = conE2 + conX2 + conU2;
conW2 = [-1 <= w2, w2 <= 1];
obj2 = O2*x2(:,2);

tau2 = sdpvar(1,1);
objP2 = tau2;
con2 = con2 + [obj2 <= tau2];

% Feasible set
conFS = conE1 + conX1 + conU1 + conW1 + conE2 + conX2 + conU2 + conW2;

%% Figure 7

% Centralized communication
dataCent= toyExampleCent(con1,con2,objP1,objP2,conW1,conW2,conFS,x1,x2,u1,u2,w1,w2,fDual,fPlot);

% Partially Nested communication
dataDecent = toyExampleDecent(con1,con2,objP1,objP2,conW1,conW2,conFS,x1,x2,u1,u2,w1,w2,fDual,fPlot);

%% Figure 8

% box set (infinity norm)
dataNinf = toyExampleRect(con1,con2,objP1,objP2,conW1,conW2,conFS,x1,x2,u1,u2,w1,w2,0,'box',fDual,fPlot,O2,xu,xl,uu,ul);

% Rombus set (one norm)
dataNone = toyExampleRect(con1,con2,objP1,objP2,conW1,conW2,conFS,x1,x2,u1,u2,w1,w2,45,'box',fDual,fPlot,O2,xu,xl,uu,ul);

% circle (2-norm)
dataNtwo = toyExampleCircle(B,C,D,O1,O2,xu,xl,uu,ul,conFS,conW1,x1,w1,0,1,fDual,fPlot);

%% Figure 9

data = cell(37,4);
for kk = 1:37
    if kk == 1 || kk == 4
        fPlot = 1;
    else
        fPlot = 0;
    end
    % rectangular sets
    dataRect = toyExampleRect(con1,con2,objP1,objP2,conW1,conW2,conFS,x1,x2,u1,u2,w1,w2,-5*(kk-1),'rect',fDual,fPlot,O2,xu,xl,uu,ul);
    data{kk,1} = dataRect;
    
    % scaled ellipsoidal sets
    dataEll = toyExampleCircle(B,C,D,O1,O2,xu,xl,uu,ul,conFS,conW1,x1,w1,5*(kk-1),1.5,fDual,fPlot);
    data{kk,2} = dataEll;
    
    % scaled ellipsoidal sets
    dataEll = toyExampleCircle(B,C,D,O1,O2,xu,xl,uu,ul,conFS,conW1,x1,w1,5*(kk-1),3,fDual,fPlot);
    data{kk,3} = dataEll;
    
    % scaled ellipsoidal sets
    dataEll = toyExampleCircle(B,C,D,O1,O2,xu,xl,uu,ul,conFS,conW1,x1,w1,5*(kk-1),10,fDual,fPlot);
    data{kk,4} = dataEll;
end

%% Figure 10
plot_objective_functions(dataCent, dataDecent, data)