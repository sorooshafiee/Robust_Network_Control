function data = toyExampleCircle(B,C,D,O1,O2,xu,xl,uu,ul,conFS,conW1,x1,w1,deg,sc,fDual,fPlot)


    % Agent 1
    W1B = [1,0,0;-1,0,0;...
        1,1,0;1,-1,0;...
        1,0,1;1,0,-1];
    h1B = [1;-1;zeros(4,1)];


    W11 = diag([1,-1,0]);
    W12 = diag([1,0,-1]);

    X1 = sdpvar(2,3);
    U1 = sdpvar(1,3);

    conE1 = [X1 == B*U1 + [zeros(2,1),C]];

    lambda = sdpvar(6,1);
    conX1 = [lambda >= 0, h1B'*lambda >= 0];
    conX1 = conX1 + [[xu-X1(1,1),-X1(1,2:3)]' == W1B'*lambda];

    lambda = sdpvar(6,1);
    conX1 = conX1 + [lambda >= 0, h1B'*lambda >= 0];
    conX1 = conX1 + [[xu-X1(2,1),-X1(2,2:3)]' == W1B'*lambda];

    lambda = sdpvar(6,1);
    conX1 = conX1 + [lambda >= 0, h1B'*lambda >= 0];
    conX1 = conX1 + [[-xl+X1(1,1),X1(1,2:3)]' == W1B'*lambda];

    lambda = sdpvar(6,1);
    conX1 = conX1 + [lambda >= 0, h1B'*lambda >= 0];
    conX1 = conX1 + [[-xl+X1(2,1),X1(2,2:3)]' == W1B'*lambda];


    lambda = sdpvar(6,1);
    conU1 = [lambda >= 0, h1B'*lambda >= 0];
    conU1 = conU1 + [[uu-U1(1,1),-U1(1,2:3)]' == W1B'*lambda];

    lambda = sdpvar(6,1);
    conU1 = conU1 + [lambda >= 0, h1B'*lambda >= 0];
    conU1 = conU1 + [[-ul+U1(1,1),U1(1,2:3)]' == W1B'*lambda];





    obj1 = O1*X1;

    tau1 = sdpvar(1,1); % worst case formulation
    objP1 = tau1;

    lambda = sdpvar(6,1);
    conObj1 = [lambda >= 0, h1B'*lambda >= 0];
    conObj1 = conObj1 + [[tau1-obj1(1,1),-obj1(1,2:3)]' == W1B'*lambda];


    % norm2 set communication
    lambda1 = sdpvar(1,2);
    alpha = sdpvar(1,1);
    beta = sdpvar(2,1);

    rd = deg2rad(deg);
    U = [cos(rd),-sin(rd);sin(rd),cos(rd)];
    Si = diag([1,sc*sc]);
    RSigma = U'*sqrt(Si)*U;
    invRSigma = U'/sqrt(Si)*U;

    alphaT = [alpha,zeros(1,2);zeros(1,2)',zeros(2,2)];
    betaT = [beta,zeros(2,2)];

    conB2 = [lambda1 >= 0];
    conB2 = conB2 + [[alphaT-lambda1(1,1)*W11-lambda1(1,2)*W12,(X1-betaT)'*RSigma';RSigma*(X1-betaT),alpha*eye(2)] >= 0];

    con1 = conE1+conX1+conU1+conObj1+conB2;


    % Agent 2
    W2 = [1,zeros(1,4);zeros(1,4)',diag([-1,-1,0,0])];
    S2 = [1,zeros(1,4);zeros(1,4)',diag([0,0,-1,-1])];



    X1T = sdpvar(2,5);
    X2 = sdpvar(2,5);
    U2 = sdpvar(1,5);

    conE2 = [X1T == [beta,zeros(2,2),alpha*invRSigma]];
    conE2 = conE2 + [X2 == B*U2 + [zeros(2,1),C,zeros(2,2)] + D*X1T];

    lambda2 = sdpvar(7,2);
    conX2 = [lambda2 >= 0];
    conX2 = conX2 + ...
        [[xu-X2(1,1),-1/2*X2(1,2:5);-1/2*X2(1,2:5)',zeros(4,4)] >= lambda2(1,1)*W2 + lambda2(1,2)*S2, ...
        [xu-X2(2,1),-1/2*X2(2,2:5);-1/2*X2(2,2:5)',zeros(4,4)] >= lambda2(2,1)*W2 + lambda2(2,2)*S2, ...
        [-xl+X2(1,1),+1/2*X2(1,2:5);+1/2*X2(1,2:5)',zeros(4,4)] >= lambda2(3,1)*W2 + lambda2(3,2)*S2, ...
        [-xl+X2(2,1),+1/2*X2(2,2:5);+1/2*X2(2,2:5)',zeros(4,4)] >= lambda2(4,1)*W2 + lambda2(4,2)*S2];


    conU2 = [[uu-U2(1,1),-1/2*U2(1,2:5);-1/2*U2(1,2:5)',zeros(4,4)] >= lambda2(5,1)*W2 + lambda2(5,2)*S2, ...
        [-ul+U2(1,1),+1/2*U2(1,2:5);+1/2*U2(1,2:5)',zeros(4,4)] >= lambda2(6,1)*W2 + lambda2(6,2)*S2];





    obj2 = O2*X2;

    tau2 = sdpvar(1,1); % worst case formulation
    objP2 = tau2;
    conObj2 = [[tau2-obj2(1,1),obj2(2:5);obj2(2:5)',zeros(4,4)] >= lambda2(7,1)*W2 + lambda2(7,2)*S2];


    con2 = conE2+conX2+conU2+conObj2;


    fprintf('\n----------- Ellipsoid rotated by %d degress and scaled by %d -----------\n',deg,sc);
    optimize(con1+con2,objP1+objP2,sdpsettings('verbose',0,'solver','sdpt3'));
    fprintf('Decentalized Problem: Robust problem objective: %f + %f = %f\n',value(objP1),...
        value(objP2),value(objP1+objP2));



    % save values
    al = value(alpha);
    bl = value(beta);

    objPrim = value(objP1+objP2);
    objPrimAg1 = value(objP1);
    objPrimAg2 = value(objP2);

    if fDual
        objDual = NaN;
    else
        objDual = NaN;
    end

    if fPlot

        % compute binding scenarios (for agent 2) by maximizing each constraint
        w2 = sdpvar(2,1);
        s2 = sdpvar(2,1);
        bindScenarios = [];

        con = [norm(w2,2) <= 1];
        con = con + [norm(s2,2) <= 1];

        obj = xu-value(X2(1,:))*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = -xl+value(X2(1,:))*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = xu-value(X2(2,:))*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = -xl+value(X2(2,:))*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = uu-value(U2(1,:))*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = -ul+value(U2(1,:))*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = value(tau2)-O2*value(X2)*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        % plot the feasible region
        x110 = value(X1(1,1));
        x111 = value(X1(1,2));
        x112 = value(X1(1,3));


        x120 = value(X1(2,1));
        x121 = value(X1(2,2));
        x122 = value(X1(2,3));

        x110D = value(X1T(1,1));
        x111D = value(X1T(1,2));
        x112D = value(X1T(1,3));
        x113D = value(X1T(1,4));
        x114D = value(X1T(1,5));


        x120D = value(X1T(2,1));
        x121D = value(X1T(2,2));
        x122D = value(X1T(2,3));
        x123D = value(X1T(2,4));
        x124D = value(X1T(2,5));

        conf = conW1 + [x1(1,2) == x110 + [x111,x112]*w1];
        conf = conf + [x1(2,2) == x120 + [x121,x122]*w1];


        conBox = [norm(RSigma*(x1(:,2)-bl),2)<= al];


        nBds = size(bindScenarios,1);
        x11bdD = zeros(nBds,1);
        x12bdD = zeros(nBds,1);
        for kk = 1:nBds
            x11bdD(kk) = [x110D,x111D,x112D,x113D,x114D]*bindScenarios(kk,:)';
            x12bdD(kk) = [x120D,x121D,x122D,x123D,x124D]*bindScenarios(kk,:)';
        end
        
        fig = figure;
        fig.Color = 'w';
        fig.PaperPositionMode = 'manual';
        fig.PaperUnits = 'inches';
        fig.Units = 'inches';
        fig.PaperPosition = [0, 0, 5, 3.09];
        fig.PaperSize = [5, 3.09];
        fig.Position = [0.5, 0.5, 4.9, 2.99];
        fig.Resize = 'off';
        fig.InvertHardcopy = 'off';

        hold on;
        plot(conFS,x1(:,2),[244, 182, 182]/255)
        plot(conBox, x1(:,2),[0.6, 0.8, 0.6]);
        plot(conf,x1(:,2),[179, 179, 249]/255)
        plot(x11bdD,x12bdD, '*k', 'MarkerSize', 8, 'LineWidth',1.25)
        grid on;

        ax = gca;
        ax.FontName = 'LaTeX';
        ax.Title.Interpreter = 'LaTeX';
        ax.XLabel.Interpreter = 'LaTeX';
        ax.YLabel.Interpreter = 'LaTeX';
        ax.TickLabelInterpreter = 'LaTeX';
        ax.Box = 'on';
        ax.XLim = [-7, 7];
        ax.XTick = [-5,0,5];
        ax.YLim = [-7, 7];
        ax.YTick = [-5,0,5];
        ax.FontSize = 14;
        ax.Title.Position = [0.5 5.5 1];

        xlabel('$x_{11}$');
        ylabel('$x_{12}$');

        axis([-7,7,-7,7]);
        if objPrimAg2 > 0
            title(sprintf('Obj. = %.2f + %.2f = %.2f\n',objPrimAg1,...
                objPrimAg2,objPrim));
        else
            title(sprintf('Obj. = %.2f - %.2f = %.2f\n',objPrimAg1,...
                objPrimAg2,objPrim));
        end
        saveas(gcf,sprintf('./Figs/FigCircle_%d_%d.pdf',round(10*sc),-deg))
    end
    data.objPrim = objPrim;
    data.objPrimAg1 = objPrimAg1;
    data.objPrimAg2 = objPrimAg2;
    data.objDual = objDual;
end