classdef RobustProblem < handle
    %----------------------------------------------------------------------
    % This class is used to robustify the problem
    %----------------------------------------------------------------------
    
    properties (SetAccess = private)
        x           % state variables
        u           % input variables ([u;uLD;uRB])
        w           % uncertain variables
        idxLD       % indices of inputs that depend on current and all past disturbances
        idxRB       % indices of inputs that are robust against disturbances
        wE          % expectation of uncertain variables
        p           % parameters
        Nhor        % horizon length
        
        originalP   % original problem
        
        robustPB    % robust problem
        
        scenarioLB  % scenario lower bound problem
    end
    
    methods
        function obj = RobustProblem(object,constr,constrW,x,u,uLD,uRB,w,wE,p)
            if nargin < 8
                error('Not all arguments were provided')
            end
            obj.originalP.object = object;    % objective
            obj.originalP.constr = constr;    % constraints
            obj.originalP.constrW = constrW;  % uncertainty set
            obj.x = x;  % states
            
            idxLD = size(u,1)+1:size(u,1)+size(uLD,1);
            idxRB = size(u,1)+size(uLD,1)+1:size(u,1)+size(uLD,1)+size(uRB,1);
            u = [u;uLD;uRB];
            
            obj.u = u;  % inputs
            obj.idxLD = idxLD; % indices of inputs that depend on current and all past disturbances
            obj.idxRB = idxRB; % indices of inputs that are robust against disturbances
            
            obj.w = w;  % uncertainty
            if isempty(wE)
                wE = zeros(size(w));
            end
            obj.wE = wE;    % expectation of uncertain variables
            
            Np = max([size(w,2),size(u,2),size(x,2)-1]);
            obj.Nhor = Np;    % horizon length
            
            if nargin < 10
                obj.p = [];  % parameters
            else
                obj.p = p;
            end
            
            % extract the problem matrices
            dataP = getProblemData(object,constr,x,u,w,Np);
            dataW = getProblemData([],constrW,[],[],w,Np);
            
            pE = w(:)'*wE(:);
            dataWE = getProblemData(pE,[],[],[],w,Np);
            
            obj.originalP.dataP = dataP;
            obj.originalP.dataW = dataW;
            obj.originalP.dataWE = dataWE;
            
            obj.robustPB = [];
            
            obj.scenarioLB = [];
        end
        
        
        
        
        function robustPB = robustify(self,flagVar,Hist)
            % function to robustify the problem
            % flagVar: 1 -> LDR depends on current uncertainty
            %          0 -> strictly causal LDR
            % Hist: -1 -> Full history
            %       >0 -> history to be considered
            
            if nargin < 2
                flagVar = 0;
                Hist = -1;
            elseif nargin < 3
                Hist = -1;
            end
            
            % dimensions
            nx = size(self.x,1);
            nu = size(self.u,1);
            nw = size(self.w,1);
            Np = self.Nhor;
            
            
            
            % Linear decision rules (LDR)
            xR = cell(nx,Np+1);
            xRv = cell(nx*(Np+1),1);  % to ease the implementation of LDR
            uR = cell(nu,Np);
            uRv = cell(nu*Np,1);
            
            
            if Hist == -1
                for kk = 1:Np+1
                    if kk == 1
                        for ii = 1:nx
                            tv = self.x(ii,1);  % x0
                            xR{ii,kk} = tv;
                            xRv{(kk-1)*nx+ii} = tv;
                        end
                    else
                        for ii = 1:nx
                            % affine policies x_t = x_{0,t} + X_t*w_hist
                            tv = sdpvar(1,1+(kk-1)*nw);
                            xR{ii,kk} = tv;
                            xRv{(kk-1)*nx+ii} = tv;
                        end
                    end
                    
                    for ii = 1:nu
                        if ~ismember(ii,self.idxRB) && ~ismember(ii,self.idxLD)
                            if kk < Np+1
                                if flagVar
                                    tv = sdpvar(1,1+kk*nw);
                                else
                                    tv = sdpvar(1,1+(kk-1)*nw);
                                end
                                uR{ii,kk} = tv;
                                uRv{(kk-1)*nu+ii} = tv;
                            end
                        elseif ismember(ii,self.idxLD)
                            if kk < Np+1
                                tv = sdpvar(1,1+kk*nw);
                                uR{ii,kk} = tv;
                                uRv{(kk-1)*nu+ii} = tv;
                            end
                        else
                            if kk < Np+1
                                tv = sdpvar(1,1);
                                uR{ii,kk} = tv;
                                uRv{(kk-1)*nu+ii} = tv;
                            end
                        end
                    end
                end
            else
                for kk = 1:Np+1
                    if kk == 1
                        for ii = 1:nx
                            tv = self.x(ii,1);
                            xR{ii,kk} = tv;
                            xRv{(kk-1)*nx+ii} = tv;
                        end
                    else
                        for ii = 1:nx
                            % affine policies x = x0 + x*w_hist
                            tv = sdpvar(1,1+(kk-1)*nw);
                            xR{ii,kk} = tv;
                            xRv{(kk-1)*nx+ii} = tv;
                        end
                    end
                    
                    for ii = 1:nu
                        if ~ismember(ii,self.idxRB) && ~ismember(ii,self.idxLD)
                            if kk < Np+1
                                if flagVar
                                    tv = [sdpvar(1,1)]; % affine term
                                    for nn = 1:kk*nw
                                        if nn <= (kk-Hist-1)*nw
                                            tv = [tv,zeros(1,1)];
                                        else
                                            tv = [tv,sdpvar(1,1)];
                                        end
                                    end
                                else
                                    tv = [sdpvar(1,1)]; % affine term
                                    for nn = 1:(kk-1)*nw
                                        if nn <= (kk-Hist-1)*nw
                                            tv = [tv,zeros(1,1)];
                                        else
                                            tv = [tv,sdpvar(1,1)];
                                        end
                                    end
                                end
                                uR{ii,kk} = tv;
                                uRv{(kk-1)*nu+ii} = tv;
                            end
                        elseif ismember(ii,self.idxLD)
                            if kk < Np+1
                                tv = sdpvar(1,1+kk*nw);
                                uR{ii,kk} = tv;
                                uRv{(kk-1)*nu+ii} = tv;
                            end
                        else
                            if kk < Np+1
                                tv = sdpvar(1,1);
                                uR{ii,kk} = tv;
                                uRv{(kk-1)*nu+ii} = tv;
                            end
                        end
                    end
                end
            end
            
            
            % variable dependency
            constrRobEq = [];
            constrRobIn = [];
            constrRobInfo = [];
            
            % create equality constraints
            Ax = self.originalP.dataP.Ax;
            Au = self.originalP.dataP.Au;
            Aw = self.originalP.dataP.Aw;
            b = self.originalP.dataP.b;
            
            tagEq = self.originalP.dataP.tagEq;
            
            neq = size(Ax,1);
            [rAx,cAx] = find(Ax);
            [rAu,cAu] = find(Au);
            
            Atw = cell(neq,Np*nw+1);    % Np*nw +1 to account for constant terms
            for key = 1:length(rAx)
                % iterate over non-empty elements of Ax (store only the
                % very necessary elements - efficient memory allocation)
                ii = rAx(key);
                kk = cAx(key);
                for jj = 1:length(xRv{kk})
                    if isempty(Atw{ii,jj})
                        Atw{ii,jj} = 0;
                    end
                    Atw{ii,jj} = Atw{ii,jj} + Ax(ii,kk)*xRv{kk}(jj);
                end
            end
            
            for key = 1:length(rAu)
                ii = rAu(key);
                kk = cAu(key);
                for jj = 1:length(uRv{kk})
                    if ~isequal(uRv{kk}(jj),0)  % case of zero for uRv
                        if isempty(Atw{ii,jj})
                            Atw{ii,jj} = 0;
                        end
                        Atw{ii,jj} = Atw{ii,jj} + Au(ii,kk)*uRv{kk}(jj);
                    end
                end
            end
            
            
            
            for ii = 1:neq
                constrRobInfo = [constrRobInfo;ConstraintInfo()];   % structure to store infromation about constraints
                constrRobInfo(end).tag = tagEq{ii};
                for jj = 1:Np*nw+1
                    if jj > 1
                        if ~isequal(Aw(ii,jj-1),0)
                            if isempty(Atw{ii,jj})
                                Atw{ii,jj} = 0;
                            end
                            Atw{ii,jj} = Atw{ii,jj} + Aw(ii,jj-1);
                        end
                        if ~isempty(Atw{ii,jj})
                            constrRobEq = [constrRobEq;Atw{ii,jj}];
                        end
                    else
                        if ~isempty(Atw{ii,jj})
                            constrRobEq = [constrRobEq;Atw{ii,jj}-b(ii)];
                            constrRobInfo(end).dualVar = length(constrRobEq);
                            constrRobInfo(end).dualVarTag = 'eq';
                        else
                            if ~isequal(b(ii),0)
                                constrRobEq = [constrRobEq;b(ii)];
                                constrRobInfo(end).dualVar = length(constrRobEq);
                                constrRobInfo(end).dualVarTag = 'eq';
                            end
                        end
                    end
                end
            end
            
            
            
            % create inequality constraints
            Fx = self.originalP.dataP.Fx;
            Fu = self.originalP.dataP.Fu;
            Fw = self.originalP.dataP.Fw;
            g = self.originalP.dataP.g;
            
            tagIn = self.originalP.dataP.tagIn;
            
            nin = size(g,1);
            [rFx,cFx] = find(Fx);
            [rFu,cFu] = find(Fu);
            
            
            
            
            % find the variables (apart from the parameters) that affect the
            % uncetainty set as a there exist constraint
            var0 = recover(setdiff(unique([getvariables(self.originalP.dataW.g),...
                getvariables(self.originalP.dataW.b)]),getvariables(self.p)));
            lvar0 = length(var0);
            
            
            % construct the uncertain matrices
            if lvar0
                % create all the relevant matrices, including the ones to
                % deal with the existence constraints
                G = self.originalP.dataW.Fw;
                [K,h,~] = getParamMatrices(self.originalP.dataW.g,var0);
                nUin = size(h,1);
                if isempty(G)
                    G = zeros(nUin,Np*nw);
                end
                if isempty(K)
                    K = zeros(nUin,lvar0);
                end
                if isempty(h)
                    h = zeros(nUin,1);
                end
                
                E = self.originalP.dataW.Aw;
                [L,f,~] = getParamMatrices(self.originalP.dataW.b,var0);
                nUeq = size(f,1);
                if isempty(E)
                    E = zeros(nUeq,Np*nw);
                end
                if isempty(L)
                    L = zeros(nUeq,lvar0);
                end
                if isempty(f)
                    f = zeros(nUeq,1);
                end
            else
                % get the uncertain matrices directly
                G = self.originalP.dataW.Fw;
                h = self.originalP.dataW.g;
                E = self.originalP.dataW.Aw;
                f = self.originalP.dataW.b;
            end
            
            
            
            Ftw = cell(nin,Np*nw+1);
            for key = 1:length(rFx)
                ii = rFx(key);
                kk = cFx(key);
                for jj = 1:length(xRv{kk})
                    if isempty(Ftw{ii,jj})
                        Ftw{ii,jj} = 0;
                    end
                    Ftw{ii,jj} = Ftw{ii,jj} + Fx(ii,kk)*xRv{kk}(jj);
                end
            end
            
            for key = 1:length(rFu)
                ii = rFu(key);
                kk = cFu(key);
                for jj = 1:length(uRv{kk})
                    if ~isequal(uRv{kk}(jj),0)  % case of zero for uRv
                        if isempty(Ftw{ii,jj})
                            Ftw{ii,jj} = 0;
                        end
                        Ftw{ii,jj} = Ftw{ii,jj} + Fu(ii,kk)*uRv{kk}(jj);
                    end
                end
            end
            
            
            % inequality contraints handling
            for ii = 1:nin
                constrRobInfo = [constrRobInfo;ConstraintInfo()];
                constrRobInfo(end).tag = tagIn{ii};
                
                Ai = [];
                idxi = [];
                for jj = 1:Np*nw+1
                    if jj > 1
                        if ~isequal(Fw(ii,jj-1),0)
                            if isempty(Ftw{ii,jj})
                                Ftw{ii,jj} = 0;
                            end
                            Ftw{ii,jj} = Ftw{ii,jj} + Fw(ii,jj-1);
                        end
                        if ~isempty(Ftw{ii,jj})
                            Ai = [Ai,Ftw{ii,jj}];
                            idxi = [idxi,jj-1];
                        end
                    else
                        bi = g(ii);
                        if ~isempty(Ftw{ii,jj})
                            bi = -Ftw{ii,jj} + bi;
                        end
                    end
                end
                
                constrRobInfo(end).BindScen_dualIdxs = idxi;
                
                if isempty(Ai)
                    if ~isnumeric(bi)
                        constrRobIn = [constrRobIn;bi];
                        constrRobInfo(end).dualVar = length(constrRobIn);
                        constrRobInfo(end).dualVarTag = 'in';
                    else
                        if 0 > bi
                            warning('Wrong formulation of problem');
                        end
                    end
                else
                    
                    if lvar0
                        % create all the dual variables to deal with
                        % existence constraints.
                        rhoVar = sdpvar(nUeq,1);
                        lambdaVar = sdpvar(nUin,1);
                        constrRobIn = [constrRobIn;lambdaVar];
                        
                        for jj = 1:length(Ai)
                            tVl = lambdaVar'*G(:,idxi(jj));
                            tVr = rhoVar'*E(:,idxi(jj));
                            tV = Ai(jj) - tVl - tVr;
                            constrRobEq = [constrRobEq;tV];
                            constrRobInfo(end).BindScen_dualVars = [constrRobInfo(end).BindScen_dualVars;length(constrRobEq)];
                            constrRobInfo(end).BindScen_dualVarsTag =[constrRobInfo(end).BindScen_dualVarsTag;'eq'];
                        end
                        
                        tVh = lambdaVar'*h;
                        tVf = rhoVar'*f;
                        tV = - tVh  - tVf + bi;
                        constrRobIn = [constrRobIn;tV];
                        constrRobInfo(end).dualVar = length(constrRobIn);
                        constrRobInfo(end).dualVarTag = 'in';
                        
                        % existence contraints
                        for jj = 1:lvar0
                            tVk = lambdaVar'*K(:,jj);
                            tVp = rhoVar'*L(:,jj);
                            tV = tVk + tVp;
                            if ~isnumeric(tV)
                                constrRobEq = [constrRobEq;tV];
                            end
                        end
                        
                        
                    else
                        % fast implementation in which only a few dual variables are created
                        idxG = [];
                        idxE = [];
                        if ~isempty(G)
                            [idxG,~] = find(G(:,idxi));
                            idxG = unique(idxG);
                        end
                        if ~isempty(E)
                            [idxE,~] = find(E(:,idxi));
                            idxE = unique(idxE);
                        end
                        
                        if ~isempty(idxG) && ~isempty(idxE)
                            
                            rhoVar = sdpvar(length(idxE),1);
                            lambdaVar = sdpvar(length(idxG),1);
                            constrRobIn = [constrRobIn;lambdaVar];
                            
                            for jj = 1:length(Ai)
                                tVl = lambdaVar'*G(idxG,idxi(jj));
                                tVr = rhoVar'*E(idxE,idxi(jj));
                                tV = Ai(jj) - tVl - tVr;
                                constrRobEq = [constrRobEq;tV];
                                constrRobInfo(end).BindScen_dualVars = [constrRobInfo(end).BindScen_dualVars;length(constrRobEq)];
                                constrRobInfo(end).BindScen_dualVarsTag =[constrRobInfo(end).BindScen_dualVarsTag;'eq'];
                            end
                            
                            tVh = lambdaVar'*h(idxG);
                            tVf = rhoVar'*f(idxE);
                            tV = - tVh - tVf + bi;
                            constrRobIn = [constrRobIn;tV];
                            constrRobInfo(end).dualVar = length(constrRobIn);
                            constrRobInfo(end).dualVarTag = 'in';
                            
                        elseif ~isempty(idxG)
                            
                            lambdaVar = sdpvar(length(idxG),1);
                            constrRobIn = [constrRobIn;lambdaVar];
                            
                            for jj = 1:length(Ai)
                                tVl = lambdaVar'*G(idxG,idxi(jj));
                                tV = Ai(jj) - tVl;
                                constrRobEq = [constrRobEq;tV];
                                constrRobInfo(end).BindScen_dualVars = [constrRobInfo(end).BindScen_dualVars;length(constrRobEq)];
                                constrRobInfo(end).BindScen_dualVarsTag =[constrRobInfo(end).BindScen_dualVarsTag;'eq'];
                            end
                            
                            tVh = lambdaVar'*h(idxG);
                            tV = - tVh + bi;
                            constrRobIn = [constrRobIn;tV];
                            constrRobInfo(end).dualVar = length(constrRobIn);
                            constrRobInfo(end).dualVarTag = 'in';
                            
                        elseif ~isempty(idxE)
                            
                            rhoVar = sdpvar(length(idxE),1);
                            for jj = 1:length(Ai)
                                tVr = rhoVar'*E(idxE,idxi(jj));
                                tV = Ai(jj) - tVr;
                                constrRobEq = [constrRobEq;tV];
                                constrRobInfo(end).BindScen_dualVars = [constrRobInfo(end).BindScen_dualVars;length(constrRobEq)];
                                constrRobInfo(end).BindScen_dualVarsTag =[constrRobInfo(end).BindScen_dualVarsTag;'eq'];
                            end
                            
                            tVf = rhoVar'*f(idxE);
                            tV =  - tVf + bi;
                            constrRobIn = [constrRobIn;tV];
                            constrRobInfo(end).dualVar = length(constrRobIn);
                            constrRobInfo(end).dualVarTag = 'in';
                            
                        else
                            error('Uncertain variables are unbounded');
                        end 
                    end
                    
                end
                
            end
            
            
            
            % create the objective
            wEt = self.originalP.dataWE.cw;
            cx = self.originalP.dataP.cx;
            cu = self.originalP.dataP.cu;
            f = self.originalP.dataP.f;
            
            ccx = find(cx);
            ccu = find(cu);
            
            ctw = cell(Np*nw+1);
            ctw(:) = {0};
            objRob = f;
            for key = 1:length(ccx)
                kk = ccx(key);
                for jj = 1:length(xRv{kk})
                    ctw{jj} = ctw{jj} + cx(kk)*xRv{kk}(jj);
                end
            end
            
            for key = 1:length(ccu)
                kk = ccu(key);
                for jj = 1:length(uRv{kk})
                    ctw{jj} = ctw{jj} + cu(kk)*uRv{kk}(jj);
                end
            end
            
            for jj = 1:Np*nw+1
                if jj > 1
                    objRob = objRob + ctw{jj}*wEt(jj-1);
                else
                    objRob = objRob + ctw{jj};
                end
            end
            
            
            % create the parametric problem
            [Feq,Peq,~] = getParamMatrices(constrRobEq,self.p);
            [Fin,Pin,paramVarN] = getParamMatrices(constrRobIn,self.p);
            
            % return matrices
            robustPB.Aeq = -Feq;
            robustPB.beq = Peq;
            robustPB.Ain = -Fin;
            robustPB.bin = Pin;
            robustPB.objRobust = objRob;
            robustPB.paramVar = paramVarN;
            robustPB.constrRobInfo = constrRobInfo;
            
            robustPB.x = xR;
            idxu = setdiff(1:nu,[self.idxLD;self.idxRB]);
            robustPB.u = uR(idxu,:);
            robustPB.uLD = uR(self.idxLD,:);
            robustPB.uRB = uR(self.idxRB,:);
            
            self.robustPB = robustPB;
        end
        
        
        
        function [conRobust,objRobust] = getRobust(self,pVal)
            % get the robust problem with respect to the parameter pVal
            
            if nargin < 2
                pVal = [];
            end
            
            
            if ~isempty(self.p)
                % replace the parameters to get the correct ordering
                paramVar = replace(self.robustPB.paramVar,self.p,pVal);
                
                % create the linear optimization problem
                conRobustEq = self.robustPB.Aeq*paramVar == self.robustPB.beq;
                conRobustIn = self.robustPB.Ain*paramVar <= self.robustPB.bin;
                
                objRobust = replace(self.robustPB.objRobust,self.p,pVal);
            else
                % create the linear optimization problem
                conRobustEq = 0 == self.robustPB.beq;
                conRobustIn = 0 <= self.robustPB.bin;
                
                objRobust = self.robustPB.objRobust;
            end
            conRobust = conRobustEq+conRobustIn;
            
            
            
            
            self.robustPB.constrEq = conRobustEq;
            self.robustPB.constrIn = conRobustIn;
            self.robustPB.constr = conRobust;
            self.robustPB.object = objRobust;
            self.robustPB.pValue = pVal;
        end
        
        
        function objectValue = solveRobust(self,solver,verbose)
            % solve the robust problem
            
            if nargin < 2
                diagnostics = optimize(self.robustPB.constr,self.robustPB.object,sdpsettings('verbose',0));
                if diagnostics.problem ~= 0
                    diagnostics = optimize(self.robustPB.constr,0,sdpsettings('verbose',0));
                    if diagnostics.problem == 0
                        error('Unbounded problem');
                    else
                        error('Infeasible problem');
                    end
                end
            else
                diagnostics = optimize(self.robustPB.constr,self.robustPB.object,sdpsettings('solver',solver,'verbose',verbose));
                if diagnostics.problem ~= 0
                    diagnostics = optimize(self.robustPB.constr,0,sdpsettings('solver',solver,'verbose',verbose));
                    if diagnostics.problem == 0
                        error('Unbounded problem');
                    else
                        error('Infeasible problem');
                    end
                end
            end
            
            
            
            
            objectValue = value(self.robustPB.object);
            fprintf('Robust problem objective: %f\n',objectValue);
            
            self.robustPB.objectValue = objectValue;
        end
        
        
        
        
        function scenario(self)
            % get the binding scenarios and formulate the lower bound
            % problem, assuming that a primal solution of the problem has
            % already been achieved.
            
            % get the dual variables
            dualEq = dual(self.robustPB.constrEq(1));
            dualIn = dual(self.robustPB.constrIn(1));
            
            % extract the scenarios
            nwF = length(self.w(:));
            bindScen = [];
            for ii = 1:length(self.robustPB.constrRobInfo)
                % constraint information
                dualVar = self.robustPB.constrRobInfo(ii).dualVar;
                dualVarTag = self.robustPB.constrRobInfo(ii).dualVarTag;
                BindScen_dualVars = self.robustPB.constrRobInfo(ii).BindScen_dualVars;
                BindScen_dualVarsTag = self.robustPB.constrRobInfo(ii).BindScen_dualVarsTag;
                BindScen_dualIdxs = self.robustPB.constrRobInfo(ii).BindScen_dualIdxs;
                
                % get the constraint dual variable
                if isequal(dualVarTag,'eq')
                    dualVar_Val = dualEq(dualVar);
                else
                    dualVar_Val = dualIn(dualVar);
                end
                
                % get the dual variables of the equality constraints
                LBindScen = length(BindScen_dualIdxs);
                BindScen_dualVars_Val = zeros(LBindScen,1);
                for kk = 1:LBindScen
                    if isequal(BindScen_dualVarsTag(kk,:),'eq')
                        BindScen_dualVars_Val(kk) = dualEq(BindScen_dualVars(kk));
                    else
                        BindScen_dualVars_Val(kk) = dualIn(BindScen_dualVars(kk));
                    end
                end
                
                % save temporary values to the constraint info
                self.robustPB.constrRobInfo(ii).dualVar_Val = dualVar_Val;
                self.robustPB.constrRobInfo(ii).BindScen_dualVars_Val = BindScen_dualVars_Val;
                
                % get the scenarios
                if dualVar_Val > 0 && LBindScen > 0
                    tVar = nan(1,nwF);
                    for kk = 1:LBindScen
                        % negative sign is needed because derivation is for
                        % >=0, while here we are using <=0. So we need -
                        % (check the paper of Chatzigiannis for details)
                        % i.e. sigma' xi >=0 for all xi in Xi then sigma in
                        % dual cone of Xi (tested that otherwise it results
                        % outside the feasible set Xi)
                        tVar(BindScen_dualIdxs(kk)) = -BindScen_dualVars_Val(kk)/dualVar_Val;
                    end
                    bindScen = [bindScen;tVar];
                end
                
                
            end
            
            % eliminate duplicated scenarios (e.g. [-1,1,NAN,NAN],[-1,1,2,NAN] -> [-1,1,2,NAN] since
            % the second scenario binds both constraints)
            Lbd = size(bindScen,1);
            idxEl = [];
            for ii = 1:Lbd
                tBd = find(~isnan(bindScen(ii,:)));
                kk = 1;
                flagEl = 0;
                while kk <= Lbd && ~flagEl
                    if kk ~= ii && ~ismember(kk,idxEl)
                        if sum(abs(bindScen(ii,tBd)-bindScen(kk,tBd))) < 1e-8
                            flagEl = 1;
                        end
                    end
                    kk = kk + 1;
                end
                if flagEl
                    idxEl = [idxEl,ii];
                end
            end
            idxNEl = setdiff([1:Lbd],idxEl);
            bindScen = bindScen(idxNEl,:);
            
            
            
            % return the binding scenarios
            self.robustPB.bindingScenarios = bindScen;
            self.scenarioLB.bindingScenarios = bindScen;
        
        end
        
        
        
        
        
        function getScenarioLB(self)
            
            wvec = self.w(:);
            nwF = length(wvec);
            bindScen = self.scenarioLB.bindingScenarios;
            
            % Extract a possible scenario:
            % If rows with nan then a feasibility problem needs to be solved to find a possible binding scenario
            Lbd = size(bindScen,1);
            for ii = 1:Lbd
                tBd = find(~isnan(bindScen(ii,:)));
                if length(tBd) < nwF
                    cons = self.originalP.constrW + [wvec(tBd) == bindScen(ii,tBd)'];
                    optimize(cons,0,sdpsettings('verbose',0));
                    bindScen(ii,:) = value(wvec);
                end
            end
            
            
            
            
            
            % Create the scenario tree - this is in order to respect the anticipativity of the problem
            % e.g. the scenarios [-1,1,2,3] and [-1,1,0,3] should have the same input variables for the
            % first two stages. (causality in our problem)
            
            nw = size(self.w,1);
            nx = size(self.x,1);
            nu = size(self.u,1);
            Np = self.Nhor;
            nrbs = size(bindScen,1);
            
            
            scnTree = [];   % stores the leafs of the tree
            scnParent = ScenarioNode();
            for bs = 1:nrbs
                pN = scnParent;
                for kk = 1:Np
                    flag = 1;
                    for ii = 1:length(pN.childVec)
                        key = pN.childVec(ii);
                        if sum(abs(key.w-bindScen(bs,(kk-1)*nw+1:kk*nw))) < 1e-8
                            pN = key;
                            flag = 0;
                            break;
                        end
                    end
                    
                    if flag
                        pN.add_child(ScenarioNode(bindScen(bs,(kk-1)*nw+1:kk*nw),nu,pN));
                        pN = pN.childVec(end);
                        if kk == Np
                            scnTree = [scnTree;pN];
                        end
                    end
                end
            end
            
            
            
            % create the optimization problem
            
            if ~isempty(self.p)
                Ax = replace(self.originalP.dataP.Ax,self.p,self.robustPB.pValue);
                Au = replace(self.originalP.dataP.Au,self.p,self.robustPB.pValue);
                Aw = replace(self.originalP.dataP.Aw,self.p,self.robustPB.pValue);
                b = replace(self.originalP.dataP.b,self.p,self.robustPB.pValue);
                
                Fx = replace(self.originalP.dataP.Fx,self.p,self.robustPB.pValue);
                Fu = replace(self.originalP.dataP.Fu,self.p,self.robustPB.pValue);
                Fw = replace(self.originalP.dataP.Fw,self.p,self.robustPB.pValue);
                g = replace(self.originalP.dataP.g,self.p,self.robustPB.pValue);
                
                cx = replace(self.originalP.dataP.cx,self.p,self.robustPB.pValue);
                cu = replace(self.originalP.dataP.cu,self.p,self.robustPB.pValue);
                f = replace(self.originalP.dataP.f,self.p,self.robustPB.pValue);
                
                x0 = replace(self.x(:,1),self.p,self.robustPB.pValue);
            else
                Ax = self.originalP.dataP.Ax;
                Au = self.originalP.dataP.Au;
                Aw = self.originalP.dataP.Aw;
                b = self.originalP.dataP.b;
                
                Fx = self.originalP.dataP.Fx;
                Fu = self.originalP.dataP.Fu;
                Fw = self.originalP.dataP.Fw;
                g = self.originalP.dataP.g;
                
                cx = self.originalP.dataP.cx;
                cu = self.originalP.dataP.cu;
                f = self.originalP.dataP.f;
                
                x0 = self.x(:,1);
            end
            
            objs = f;
            cons = [];
            
            nrbs = size(scnTree,1);
            for ii = 1:nrbs
                
                xtemp = [x0;sdpvar(nx*Np,1)];
                utemp = [];
                wtemp = [];
                
                % reverse backwards the tree
                pN = scnTree(ii);
                for kk = 1:Np
                    utemp = [pN.u;utemp];
                    wtemp = [pN.w;wtemp];
                    pN = pN.parentNode;
                end
                
                cons = cons + [Ax*xtemp + Au*utemp + Aw*wtemp == b];
                cons = cons + [Fx*xtemp + Fu*utemp + Fw*wtemp <= g];
                
                objs = objs + 1/nrbs*cx*xtemp + 1/nrbs*cu*utemp;
            end
            
            % return the binding scenarios and constraints objective
            self.scenarioLB.constr = cons;
            self.scenarioLB.object = objs;
        end
        
        
        
        
        
        function objectValue = solveScenarioLB(self,solver,verbose)
            % solve the scenario problem
            
            if nargin < 2
                optimize(self.scenarioLB.constr,self.scenarioLB.object,sdpsettings('verbose',0));
            else
                optimize(self.scenarioLB.constr,self.scenarioLB.object,sdpsettings('solver',solver,'verbose',verbose));
            end
            
            
            
            
            objectValue = value(self.scenarioLB.object);
            fprintf('Scenario lower bound objective: %f\n',objectValue);
            
            self.scenarioLB.objectValue = objectValue;
        end
        
        
        
        
        function dualCon = getDualConstraint(self,conTag)
            
            constrRobInfo = self.robustPB.constrRobInfo; 
            Lcon = length(constrRobInfo);
            
            dualConT = [];
            for kk = 1:Lcon
                if isequal(constrRobInfo(kk).tag,conTag)
                    dualConT = [dualConT;constrRobInfo(kk).dualVar_Val];
                end
            end
            
            % return value
            dualCon = dualConT;
        end
        
        
    end
    
    
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%       FUNCTIONS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = getProblemData(object,constr,x,u,w,Np)
% Returns the problem's matrices
%
% minimize cx(p) x + cu(p) u + cw(p) w + f(p)
% subj. to Ax(p) x + Au(p) u + Aw(p) w == b(p)
%          Fx(p) x + Fu(p) u + Fw(p) w <= g(p)


% validity check
if isempty(object)
    object = 0;
end
if isempty(constr)
    constr = [];
end

% get dimensions
nx = size(x,1);
nu = size(u,1);
nw = size(w,1);


% problem variables
vars = [x(:);u(:);w(:)];

% get the equality and inequality constaints
F_eq = extractConstraints(constr,'equality');
LFeq = length(F_eq);
if LFeq
    F_in = constr - F_eq;
else
    F_in = constr;
end
LFin = length(F_in);


% You have to split the equalities and inequalies because the sdpvar
% command works for con >= 0 and do not handle the equalities properly
% (please be carefull, the way it is used here is the proper).

FeqVar = [];
tagEq = [];
FinVar = [];
tagIn = [];
if LFeq
    for kk = 1:LFeq
        tVar = sdpvar(F_eq(kk));
        FeqVar = [FeqVar;tVar];
        for nn = 1:length(tVar)
            tagEq = [tagEq;{tag(F_eq(kk))}];
        end
    end
end
if LFin
    for kk = 1:LFin
        tVar = sdpvar(F_in(kk));
        FinVar = [FinVar;tVar];
        for nn = 1:length(tVar)
            tagIn = [tagIn;{tag(F_in(kk))}];
        end
    end
end

FobjVar = object;


[c,f,~] = getParamMatrices(FobjVar,vars);
[A,b,~] = getParamMatrices(FeqVar,vars);
[F,g,~] = getParamMatrices(FinVar,vars);

% assign the correct signs, the sdpvar command gives the expression >= 0...
A = -A;
F = -F;

% extract the system matrices
idxtot = getvariables(vars);

cx = [];
cu = [];
cw = [];
Ax = [];
Au = [];
Aw = [];
Fx = [];
Fu = [];
Fw = [];
for ii = 1:Np+1
    
    if nx
        % get the variables in the same order as given by the user
        idx = [];
        for nn = 1:nx
            idx = [idx,getvariables(x(nn,ii))];
        end
        
        % objective
        for jj = 1:length(idx)
            cx = [cx,c(:,idx(jj)==idxtot)];
        end
        
        % equality constraints
        if ~isempty(A)
            for jj = 1:length(idx)
                Ax = [Ax,A(:,idx(jj)==idxtot)];
            end
        end
        
        % inequality constraints
        if ~isempty(F)
            for jj = 1:length(idx)
                Fx = [Fx,F(:,idx(jj)==idxtot)];
            end
        end
    end
    
    
    if nu && ii < Np+1
        idx = [];
        for nn = 1:nu
            idx = [idx,getvariables(u(nn,ii))];
        end
        
        for jj = 1:length(idx)
            cu = [cu,c(:,idx(jj)==idxtot)];
        end
        
        if ~isempty(A)
            for jj = 1:length(idx)
                Au = [Au,A(:,idx(jj)==idxtot)];
            end
        end
        
        if ~isempty(F)
            for jj = 1:length(idx)
                Fu = [Fu,F(:,idx(jj)==idxtot)];
            end
        end
    end
    
    
    if nw && ii < Np+1
        idx = [];
        for nn = 1:nw
            idx = [idx,getvariables(w(nn,ii))];
        end
        for jj = 1:length(idx)
            cw = [cw,c(:,idx(jj)==idxtot)];
        end
        
        if ~isempty(A)
            for jj = 1:length(idx)
                Aw = [Aw,A(:,idx(jj)==idxtot)];
            end
        end
        
        if ~isempty(F)
            for jj = 1:length(idx)
                Fw = [Fw,F(:,idx(jj)==idxtot)];
            end
        end
    end
end







% return values
data.cx = cx;
data.cu = cu;
data.cw = cw;
data.f = f;
data.Ax = Ax;
data.Au = Au;
data.Aw = Aw;
data.b = b;
data.Fx = Fx;
data.Fu = Fu;
data.Fw = Fw;
data.g = g;

data.tagEq = tagEq;
data.tagIn = tagIn;
end



function [F,G,details] = getParamMatrices(P,xReal)
% Purpose: given an expression P and the parameters x that it depends on
%          return the matrices P = Fx + G (parametric on all the other
%          variables that may be present in your problem)
%          For fast implementation provide depRows (rows that depend on
%          xReal)
% Author:  Darivianakis Georgios
%
% email:   gdarivia@control.ee.ethz.ch

if nargin < 2
    xReal = recover(getvariables(P));
end

realVar = getvariables(xReal);
rVlen = length(realVar);
% Make sure they are sorted
x = recover(realVar);
Qtemp = [];
htemp = [];
for ii = 1:length(P)
    if sum(ismember(depends(P(ii)),realVar)) == 0
        Qtemp = [Qtemp;zeros(1,rVlen)];
        htemp = [htemp;P(ii)];
    else
        
        % much faster code...
        [Q,c,f,allvar,~] = quaddecompGD(P(ii),x);
        if ~isreal(Q) % Numerical noise common on imaginary parts
            Qr = real(Q);
            Qi = imag(Q);
            Qr(abs(Qr)<1e-10) = 0;
            Qi(abs(Qi)<1e-10) = 0;
            cr = real(c);
            ci = imag(c);
            cr(abs(cr)<1e-10) = 0;
            ci(abs(ci)<1e-10) = 0;
            Q = Qr + sqrt(-1)*Qi;
            c = cr + sqrt(-1)*ci;
        end
        
        
        used_variables = getvariables(allvar);
        
        notparameters = find(ismember(used_variables,getvariables(x)));%parameters);
        parameters = setdiff(1:length(used_variables),notparameters);
        if ~isempty(parameters)
            y = recover(used_variables(parameters));
            
            % parametric case (nonlinear Pparam)
            ht = f + c(parameters)'*y + y'*Q(parameters,parameters)*y;
            c(parameters) = [];
            Q2 = Q(notparameters,parameters);
            
            c = c + 2*Q2*y;
        else
            ht = f;
        end
        
        %% create the matrices
        Qtemp = [Qtemp;c'];
        htemp = [htemp;ht];
    end
end

F = Qtemp;
G = htemp;
details = x;

end



function [Q,c,f,x,info] = quaddecompGD(p,z)
%QUADDECOMP Internal function to decompose quadratic expression

[n,m]=size(p);
info = 0;

% Is it a scalar
if (n*m==1)
    % Involved in polynomial expression
    [mt,variabletype] = yalmip('monomtable');
    x_lin = getvariables(p);
    x_var = find(any(mt(x_lin,:),1));
    if nargin==2
        x_var = union(x_var,depends(z));
    end
    x = recover(x_var);
    if all(variabletype(x_lin) ==0)% is(p,'linear')
        n = length(x);
        Q = spalloc(n,n,0);
        fc = getbase(p);
        f = fc(1);
        if nargin==2
            vars = getvariables(p);
            c = zeros(length(x),1);
            for i = 1:length(vars)
                c(vars(i)==x_var) = fc(1+i);
            end
        else
            c = fc(2:end);c=c(:);
        end
        return
    end
    variabletype = variabletype(x_lin);
    if all(variabletype<=2)
        
        base = getbase(p);
        if nnz(base(1))==0
            f = 0;
            base = base(2:end);
        else
            f = base(1);
            base = base(2:end);
        end
        mt = mt(x_lin,x_var);
        quads   = find (variabletype == 2);
        bilins  = find (variabletype == 1);
        linears  = find (variabletype == 0);
        [varsC,~,~] = find(mt(linears,:)');
        [varsQ,~,~] = find(mt(quads,:)');
        [varsB,~,~] = find(mt(bilins,:)');
        if isempty(varsQ)
            varsQ = [];
        end
        if isempty(varsB)
            varsB = [];
        end
        if isempty(varsC)
            varsC = [];
        end
        c = sparse(varsC,1,base(linears),length(x),1);
        ii = [varsQ ; varsB(1:2:end) ; varsB(2:2:end)];
        jj = [varsQ ; varsB(2:2:end) ; varsB(1:2:end)];
        kk = [base(quads)  base(bilins)/2  base(bilins)/2];
        Q = sparse(ii,jj,kk,length(x),length(x));
    else
        if nargout==5
            info = 1;
            Q = [];
            c = [];
            f = [];
            x = [];
        else
            error('Function is not quadratic');
        end
    end
    
else
    if nargout==5
        info = 1;
        Q = [];
        c = [];
        f = [];
        x = [];
    else
        error('Function is not scalar');
    end
end

end