% Purpose: Derive the robust counterpart of the optimization problem that
%          is defined by an objective function and a constraint set. 
% Author:  Darivianakis Georgios
% email:   gdarivia@control.ee.ethz.ch
% Date:    22/2/2015


function [optRobust,realVarRobust] = getRobust(con,obj,inputVar,idxNDep,idxLDep,stateVar,...
    uncertainVar,uncertainMean,uncertainUpper,uncertainLower,uncertainSlack,...
    uncertainflag,paramVar)
% [conRobust,objRobust] = getRobust(con,obj,uncertainPar,uncertainSet)
% con: constaint set of the original problem
% obj: objective of the original problem
% realVar: list of sdpvar variables that do not depend on uncertainty. 
% uncertainVar: list of the uncertain sdpvar variables
% uncertainMean: Mean value of the uncertainties (used to compute the
%                expectation of the objective function)
% uncertainSet: set at which the uncertainty belongs
% uncertainParam: slackness set for the linear decision rules that denotes
%                 the number of previous uncertainty that it depends on.
%                 For example, if we have 0: Robust, 1: w_{k-1}
% flag: = 1 : decision rule does not depend on current disturbance (not
%             realized before the control is evaluated)
%       = 0 : decision rule depends on current disturbance (realized before
%             the control action is evaluated)
% conRobust: robustified constrained set
% objRobust: robustified objective function


% get a row vector for the uncertain and real parameters

realVar = [inputVar;stateVar];
nRunc = size(uncertainVar,1);
nCunc = size(uncertainVar,2);
nRreal = size(realVar,1);
nCreal = size(realVar,2);
nRrealInput = size(inputVar,1);
nCrealInput = size(inputVar,2);
nRrealState = size(stateVar,1);
nCrealState = size(stateVar,2);
% error check
if ~isequal(nCreal,nCunc)
    error('Real and uncertain variables have not the same length');
end

%% get the matrices


% get the equality and inequality constaints

F_eq = extractConstraints(con,'equality');
F_in = con - F_eq;

% You have to split the equalities and inequalies because the sdpvar
% command works for con >= 0 and do not handle the equalities properly
% (please be carefull, the way it is used here is the proper).
FeqVar = sdpvar(F_eq);
FinVar = sdpvar(F_in);

[Few,Pe,~] = getParamMatrices(FeqVar,uncertainVar);
[Fcw,Pc,~] = getParamMatrices(FinVar,uncertainVar);

% assign the correct signs, the sdpvar command gives the expression >= 0...
Few = -Few;
Fcw = -Fcw;
Pobj = obj;

%% Assign the variables (Linear Decision Rules)

% For the same disturbances assign the same control variable.
uncertainVarReal = [];
for ii = 1:nCunc
    uncertainVarReal = [uncertainVarReal,recover(getvariables(uncertainVar(:,ii)))];
end
nRuncReal = size(uncertainVarReal,1);

% create the decision rules for the real variables (Therefore, extended
% matrices depend only on the size of the real variables)
realVarExt = cell(nRreal,nCreal);
realVarIdent = cell(nRreal,nCreal);
identS = 1110;  % starting value for identifiers
ListVarIdent = [];
uncertainVarExt = cell(nRreal,nCreal);

if uncertainflag == 0          % decision rule is not affected by the current uncertainty (uncertainty is realized after the control action)
    for jj = 1:nRrealInput
        % idxNDep : variables that are robust
        % idxLDep : variables that are being affected by current uncertainty
        if sum(sum(ismember(idxNDep,jj))) == 0 && sum(sum(ismember(idxLDep,jj))) == 0
            for kk = 1:nCrealInput
                cntVar = min([(uncertainSlack)*nRuncReal+1,(kk-1)*nRuncReal+1]);
                realVarExt{jj,kk} = sdpvar(cntVar,1);
                realVarIdent{jj,kk} = [identS:identS+cntVar-1]';
                identS = identS+cntVar;
                ListVarIdent = [ListVarIdent;[realVarIdent{jj,kk},realVarExt{jj,kk}]];
                uncertainVarExt{jj,kk} = 1;
                slackStrt = max([kk-uncertainSlack,1]);
                for ll = slackStrt:kk-1
                    uncertainVarExt{jj,kk} = [uncertainVarExt{jj,kk};uncertainVarReal(:,ll)];
                end
            end
        elseif sum(sum(ismember(idxLDep,jj))) == 1
            for kk = 1:nCrealInput
                cntVar = min([(uncertainSlack+1)*nRuncReal+1,(kk)*nRuncReal+1]);  % change of index to account and the current disturbance
                realVarExt{jj,kk} = sdpvar(cntVar,1);
                realVarIdent{jj,kk} = [identS:identS+cntVar-1]';
                identS = identS+cntVar;
                ListVarIdent = [ListVarIdent;[realVarIdent{jj,kk},realVarExt{jj,kk}]];
                uncertainVarExt{jj,kk} = 1;
                slackStrt = max([kk-uncertainSlack,1]);
                for ll = slackStrt:kk   % change of index to account and the current disturbance
                    uncertainVarExt{jj,kk} = [uncertainVarExt{jj,kk};uncertainVarReal(:,ll)];
                end
            end
        else
            % variables that do not depend on previous uncertainties
            % (Robust version)
            for kk = 1:nCrealInput
                cntVar = 1;
                realVarExt{jj,kk} = inputVar(jj,kk);
                realVarIdent{jj,kk} = [identS:identS+cntVar-1]';
                identS = identS+cntVar;
                ListVarIdent = [ListVarIdent;[realVarIdent{jj,kk},realVarExt{jj,kk}]];
                uncertainVarExt{jj,kk} = 1;
            end
        end
    end
    % the states should also depend on the current and all the previous
    % uncertainties
    for jj = nRrealInput+1:nRrealInput+nRrealState
        for kk = 1:nCrealState
            % change of index to account and the current disturbance and
            % one more to get the value of x(k-1)
            cntVar = (kk)*nRuncReal+1;  
            realVarExt{jj,kk} = sdpvar(cntVar,1);
            realVarIdent{jj,kk} = [identS:identS+cntVar-1]';
            identS = identS+cntVar;
            ListVarIdent = [ListVarIdent;[realVarIdent{jj,kk},realVarExt{jj,kk}]];
            uncertainVarExt{jj,kk} = 1;
            slackStrt = 1; % account one before to be consistent on constraints
            for ll = slackStrt:kk   % change of index to account and the current disturbance
                uncertainVarExt{jj,kk} = [uncertainVarExt{jj,kk};uncertainVarReal(:,ll)];
            end
        end
    end
   
else  % decision rule is being affected by the current uncertainty (uncertainty is realized before the control action)
    for jj = 1:nRrealInput
        if sum(sum(ismember(idxNDep,jj))) == 0
            for kk = 1:nCrealInput
                cntVar = min([(uncertainSlack+1)*nRuncReal+1,(kk)*nRuncReal+1]);  % change of index to account and the current disturbance
                realVarExt{jj,kk} = sdpvar(cntVar,1);
                realVarIdent{jj,kk} = [identS:identS+cntVar-1]';
                identS = identS+cntVar;
                ListVarIdent = [ListVarIdent;[realVarIdent{jj,kk},realVarExt{jj,kk}]];
                uncertainVarExt{jj,kk} = 1;
                slackStrt = max([kk-uncertainSlack,1]);
                for ll = slackStrt:kk   % change of index to account and the current disturbance
                    uncertainVarExt{jj,kk} = [uncertainVarExt{jj,kk};uncertainVarReal(:,ll)];
                end
            end
        else
            % variables that do not depend on previous uncertainties
            % (Robust version)
            for kk = 1:nCrealInput
                cntVar = 1;
                realVarExt{jj,kk} = inputVar(jj,kk);
                realVarIdent{jj,kk} = [identS:identS+cntVar-1]';
                identS = identS+cntVar;
                ListVarIdent = [ListVarIdent;[realVarIdent{jj,kk},realVarExt{jj,kk}]];
                uncertainVarExt{jj,kk} = 1;
            end
        end
    end
    for jj = nRrealInput+1:nRrealInput+nRrealState
        for kk = 1:nCrealState
            cntVar = kk*nRuncReal+1;  % change of index to account and the current disturbance
            realVarExt{jj,kk} = sdpvar(cntVar,1);
            realVarIdent{jj,kk} = [identS:identS+cntVar-1]';
            identS = identS+cntVar;
            ListVarIdent = [ListVarIdent;[realVarIdent{jj,kk},realVarExt{jj,kk}]];
            uncertainVarExt{jj,kk} = 1;
            slackStrt = 1;
            for ll = slackStrt:kk   % change of index to account and the current disturbance
                uncertainVarExt{jj,kk} = [uncertainVarExt{jj,kk};uncertainVarReal(:,ll)];
            end
        end
    end
end

%% Create the decision rules matrices

% Decision rules assignment with a constraint set
conLocal = [];
for jj = 1:nRreal
    for kk = 1:nCreal
        conLocal = conLocal + [realVar(jj,kk) == realVarIdent{jj,kk}'*uncertainVarExt{jj,kk}];
    end
end

conLocal = sdpvar(conLocal);

[FeL,PeL] = getParamMatrices(conLocal,realVar);

FeL = -FeL;
[FewL,GewL] = getParamMatrices(PeL,uncertainVar);


% assign the correct control variables to the positions of the identifiers
nRw = size(FewL,1);
nCw = size(FewL,2);
Fwfull = [];
Gwfull = [];
for ii = 1:nRw
    if GewL(ii)>0
        idxr = GewL(ii) - ListVarIdent(1,1)+1;
        Gwfull = [Gwfull;ListVarIdent(idxr,2)];
    else
        idxr = -GewL(ii) - ListVarIdent(1,1)+1;
        Gwfull = [Gwfull;-ListVarIdent(idxr,2)];
    end
    
    Fwtemp = [];
    for jj = 1:nCw
        if ~isequal(FewL(ii,jj),0)
            if FewL(ii,jj) > 0
                idxr = FewL(ii,jj) - ListVarIdent(1,1)+1;
                Fwtemp = [Fwtemp,ListVarIdent(idxr,2)];
            else
                idxr = -FewL(ii,jj) - ListVarIdent(1,1)+1;
                Fwtemp = [Fwtemp,-ListVarIdent(idxr,2)];
            end
        else
            Fwtemp = [Fwtemp,0];
        end
    end
    Fwfull = [Fwfull;Fwtemp];
end

Fufull = FeL;
% Fufull is invertible (make Fufull the identity matrix)
Fufull = full(Fufull);
if rank(Fufull)<size(Fufull,2)
    error('Rank problem: Matrix can not be inversed');
end

% get the inverse matrix
Fuinv = inv(Fufull);
Fw = Fuinv*Fwfull;
Gw = Fuinv*Gwfull;

%% create the correct matrices
[Fc,Gc,~] = getParamMatrices(Pc,realVar);
[Fe,Ge,~] = getParamMatrices(Pe,realVar);
[Fobj,Gobj,~] = getParamMatrices(Pobj,realVar);

%% get the final optimization problem matrices
ctot = Fobj*Fw;
htot = Fobj*Gw + Gobj;

Aetot = Few - Fe*Fw;
Betot = Fe*Gw + Ge;
Atot = Fcw - Fc*Fw;
Btot = Fc*Gw + Gc;

%% get the mean value of the disturbances as a vector of the disturbances

Pmean = 0;
PSumDiv = 0;
for jj = 1:nRunc
    for kk = 1:nCunc
        Pmean = Pmean + uncertainVar(jj,kk)*uncertainMean(jj,kk);
        PSumDiv = PSumDiv + uncertainVar(jj,kk);
    end
end
[vptMean,~,~] = getParamMatrices(Pmean,uncertainVar);
[vptSumDiv,~] = getParamMatrices(PSumDiv,uncertainVar);
if ~isempty(find(vptSumDiv==0,1))
    error('This formulation is wrong');
end
vptMean = vptMean'./vptSumDiv'; % vector of the same size as the disturbances

% get the upper bound as a vector of the disturbances

Pupper = 0;
for jj = 1:nRunc
    for kk = 1:nCunc
        Pupper = Pupper + uncertainVar(jj,kk)*uncertainUpper(jj,kk);
    end
end
[vptUpper,~,~] = getParamMatrices(Pupper,uncertainVar);
vptUpper = vptUpper'./vptSumDiv'; % vector of the same size as the disturbances

% get the lower bound as a vector of the disturbances

Plower = 0;
for jj = 1:nRunc
    for kk = 1:nCunc
        Plower = Plower + uncertainVar(jj,kk)*uncertainLower(jj,kk);
    end
end
[vptLower,~,~] = getParamMatrices(Plower,uncertainVar);
vptLower = vptLower'./vptSumDiv'; % vector of the same size as the disturbances


%% create the disturbance set matrices (Gw<=F) and mean value
Gtot = [];
Ftot = [];
for ii = 1:length(vptMean)
    Gtt = [1;-1];
    Ftt = [vptUpper(ii);-vptLower(ii)];
    Gtot = blkdiag(Gtot,Gtt);
    Ftot = [Ftot;Ftt];
end


%% create the finite dimensional robust optimization problem
vpt = vptMean;
objRobustTemp = ctot*vpt + htot;
FinRobust = [];
FeqRobust = [];

cntEq = 0;
cntIn = 0;
EqVar = [];
InVar = [];
paramVarDeps = depends(paramVar);

% equality constraints
idxEqA = [];
cnt = 0;
for jj = 1:size(Aetot,2)
    for ii = 1:size(Aetot,1)
        cnt = cnt + 1;
        if ~isnumeric(Aetot(ii,jj))
            idxEqA = [idxEqA;cnt];
        end
    end
end

for ii = 1:length(idxEqA)
    cntEq = cntEq + 1;
    FeqRobust = [FeqRobust;Aetot(idxEqA(ii))];
    if sum(ismember(depends(Aetot(idxEqA(ii))),paramVarDeps))>0
        EqVar = [EqVar;cntEq];
    end
end

idxEqB = [];
cnt = 0;
for jj = 1:size(Betot,2)
    for ii = 1:size(Betot,1)
        cnt = cnt + 1;
        if ~isnumeric(Betot(ii,jj))
            idxEqB = [idxEqB;cnt];
        end
    end
end
for ii = 1:length(idxEqB)
    cntEq = cntEq + 1;
    FeqRobust = [FeqRobust;Betot(idxEqB(ii))];
    if sum(ismember(depends(Betot(idxEqB(ii))),paramVarDeps))>0
        EqVar = [EqVar;cntEq];
    end
end


% inequality constraints
lambdaFull = sdpvar(size(Gtot,1),size(Atot,1),'full');

for ii = 1:size(Atot,1)
    if isequal(Atot(ii,:),zeros(size(Atot(ii,:))))
        if ~isnumeric(Btot(ii))
            cntIn = cntIn + 1;
            FinRobust = [FinRobust;Btot(ii)];
            if sum(ismember(depends(Btot(ii)),paramVarDeps))>0
                InVar = [InVar;cntIn];
            end
        else
            if 0 > Btot(ii)
                warning('Wrong Formulation of Problem\n')
            end
        end
    else
        % check the entries of Atot that are non-zero and create the
        % appropiate dual variables
        idx = [];
        for nn = 1:length(Atot(ii,:))
            if ~isequal(Atot(ii,nn),0)
                idx = [idx;nn];
            end
        end
        idxG = [];

        for nn = 1:length(idx)
            idxt = find(Gtot(:,idx(nn))~=0);
            idxG = union(idxG,idxt);
        end
        Gtemp = Gtot(idxG,:);
        Ftemp = Ftot(idxG);
        
        lambda = lambdaFull(1:length(idxG),ii); % dual variables for the inequality constraints of the constraint set
        
        cntIn = cntIn + 1;
        FinRobust = [FinRobust;Btot(ii)-lambda'*Ftemp;lambda];
        if sum(ismember(depends(Btot(ii)-lambda'*Ftemp),paramVarDeps))>0
            InVar = [InVar;cntIn];
        end
        cntIn = cntIn + length(idxG);   % because of lambda variables above
        
        for nn = 1:size(Atot(ii,:)',1)
            cntEq = cntEq + 1;
            FeqRobust = [FeqRobust;Atot(ii,nn)'- Gtemp(:,nn)'*lambda];
            if sum(ismember(depends(Atot(ii,nn)'- Gtemp(:,nn)'*lambda),paramVarDeps))>0
                EqVar = [EqVar;cntEq];
            end
        end
    end
end

%% return values
objRobust = objRobustTemp;

FeqVar = (FeqRobust); % do that to remove the unnecessary zeros
FinVar = (FinRobust);

[Feq,Peq,~] = getParamMatrices(FeqVar,paramVar);
[Fin,Pin,paramVarN] = getParamMatrices(FinVar,paramVar);

% return matrices
optRobust.Aeq = -Feq;
optRobust.beq = Peq;
optRobust.Ain = -Fin;
optRobust.bin = Pin;
optRobust.objRobust = objRobust;
optRobust.paramVar = paramVarN;

realVarRobust = realVarExt;
end