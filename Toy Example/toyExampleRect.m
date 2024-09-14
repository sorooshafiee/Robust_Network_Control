function data = toyExampleRect(con1,con2,objP1,objP2,conW1,conW2,conFS,x1,x2,u1,u2,w1,w2,deg,type,fDual,fPlot,O2,xu,xl,uu,ul)

    % auxiliary variable
    Np = size(u1,2);


    s2 = sdpvar(2,Np,'full'); % primitive variable
    if strcmp(type,'box')
        a1t = sdpvar(1,Np,'full'); % scaling
        a1 = [a1t;a1t];
    else
        a1 = sdpvar(2,Np,'full'); % scaling
    end
    b1 = sdpvar(2,Np,'full'); % translation

    rd = deg2rad(deg);

    % Agent 1
    conP1 = con1;
    conS1 = conW1;
    for kk = 1:Np
        conP1 = conP1 + [cos(rd)*x1(1,kk+1)-cos(rd)*b1(1,kk)+sin(rd)*x1(2,kk+1)-sin(rd)*b1(2,kk) <= a1(1,kk),...
            sin(rd)*x1(1,kk+1)-sin(rd)*b1(1,kk)-cos(rd)*x1(2,kk+1)+cos(rd)*b1(2,kk) <= a1(2,kk),...
            -sin(rd)*x1(1,kk+1)+sin(rd)*b1(1,kk)+cos(rd)*x1(2,kk+1)-cos(rd)*b1(2,kk) <= a1(2,kk),...
            -cos(rd)*x1(1,kk+1)+cos(rd)*b1(1,kk)-sin(rd)*x1(2,kk+1)+sin(rd)*b1(2,kk) <= a1(1,kk)];

    end


    robustP1 = RobustProblem(objP1,conP1,conS1,x1,u1,[],[],w1,[]);
    robustP1.robustify(1);
    robustP1.getRobust;

    % Agent 2

    conP2 = con2;
    conS2 = conW2;
    for kk = 1:Np
        conP2 = conP2 + [x1(:,kk+1) == [cos(rd)*a1(1,kk),-sin(rd)*a1(2,kk);...
            sin(rd)*a1(1,kk),cos(rd)*a1(2,kk)]*s2(:,kk)+b1(:,kk)];

        conS2 = conS2 + [-1 <= s2(:,kk), s2(:,kk) <= 1];
    end


    robustP2 = RobustProblem(objP2,conP2,conS2,x2,u2,x1(:,2:Np+1),[],[w2;s2],[]);
    robustP2.robustify(1);
    robustP2.getRobust;

    objRob = robustP1.robustPB.object + robustP2.robustPB.object;
    conRob = robustP1.robustPB.constr + robustP2.robustPB.constr;

    optimize(conRob,objRob,sdpsettings('verbose',0));
    fprintf('\n********************** Rotated %s set by %d degress **********************\n',type,deg);
    fprintf('Decentalized Problem: Robust problem objective: %f + %f = %f\n',value(robustP1.robustPB.object),...
        value(robustP2.robustPB.object),value(objRob));



    % save values
    al = value(a1(:,1));
    bl = value(b1(:,1));

    objPrim = value(objRob);
    objPrimAg1 = value(robustP1.robustPB.object);
    objPrimAg2 = value(robustP2.robustPB.object);

    if fDual
        robustP1.scenario;
        robustP2.scenario;

        robustP1.getScenarioLB;
        robustP2.getScenarioLB;

        objScen = robustP1.scenarioLB.object + robustP2.scenarioLB.object;
        conScen = robustP1.scenarioLB.constr + robustP2.scenarioLB.constr;

        optimize(conScen,objScen,sdpsettings('verbose',0));
        fprintf('Decentalized Problem: Scenario lower bound objective: %f\n\n',value(objScen));

        objDual = value(objScen);
    else
        objDual = NaN;
    end

    % plots
    if fPlot

        % compute binding scenarios (for agent 2) by maximizing each constraint
        w2 = sdpvar(2,1);
        s2 = sdpvar(2,1);
        bindScenarios = [];

        con = [-1 <= w2 <= 1];
        con = con + [-1 <= s2 <= 1];

        obj = xu-value(robustP2.robustPB.x{1,2})*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = -xl+value(robustP2.robustPB.x{1,2})*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = xu-value(robustP2.robustPB.x{2,2})*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = -xl+value(robustP2.robustPB.x{2,2})*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = uu-value(robustP2.robustPB.u{1,1})*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = -ul+value(robustP2.robustPB.u{1,1})*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        obj = value(robustP2.robustPB.object)...
            -O2*value([robustP2.robustPB.x{1,2};robustP2.robustPB.x{2,2}])*[1;w2;s2];
        optimize(con,obj,sdpsettings('verbose',0));
        if abs(value(obj))<=1e-6
            bindScenarios = [bindScenarios;value([1;w2;s2]')];
        end

        % plot the feasible region
        x110 = value(robustP1.robustPB.x{1,2}(1));
        x111 = value(robustP1.robustPB.x{1,2}(2));
        x112 = value(robustP1.robustPB.x{1,2}(3));

        x110D = value(robustP2.robustPB.uLD{1,1}(1));
        x111D = value(robustP2.robustPB.uLD{1,1}(2));
        x112D = value(robustP2.robustPB.uLD{1,1}(3));
        x113D = value(robustP2.robustPB.uLD{1,1}(4));
        x114D = value(robustP2.robustPB.uLD{1,1}(5));

        x120 = value(robustP1.robustPB.x{2,2}(1));
        x121 = value(robustP1.robustPB.x{2,2}(2));
        x122 = value(robustP1.robustPB.x{2,2}(3));

        x120D = value(robustP2.robustPB.uLD{2,1}(1));
        x121D = value(robustP2.robustPB.uLD{2,1}(2));
        x122D = value(robustP2.robustPB.uLD{2,1}(3));
        x123D = value(robustP2.robustPB.uLD{2,1}(4));
        x124D = value(robustP2.robustPB.uLD{2,1}(5));

        conf = conW1 + conW2 + [x1(1,2) == x110 + [x111,x112]*w1];
        conf = conf + [x1(2,2) == x120 + [x121,x122]*w1];

        nBds = size(bindScenarios,1);
        x11bdD = zeros(nBds,1);
        x12bdD = zeros(nBds,1);
        for kk = 1:nBds
            x11bdD(kk) = [x110D,x111D,x112D,x113D,x114D]*bindScenarios(kk,:)';
            x12bdD(kk) = [x120D,x121D,x122D,x123D,x124D]*bindScenarios(kk,:)';
        end


        conBox = [cos(rd)*(x1(1,2)-bl(1,1)) + sin(rd)*(x1(2,2)-bl(2,1)) <= al(1,1),...
            sin(rd)*(x1(1,2)-bl(1,1)) - cos(rd)*(x1(2,2)-bl(2,1)) <= al(2,1),...
            -sin(rd)*(x1(1,2)-bl(1,1)) + cos(rd)*(x1(2,2)-bl(2,1)) <= al(2,1),...
            -cos(rd)*(x1(1,2)-bl(1,1)) - sin(rd)*(x1(2,2)-bl(2,1)) <= al(1,1)];

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
   

        plot(conFS,x1(:,2),[244, 182, 182]/255);
        hold on;
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
        if value(objP2) > 0
            title(sprintf('Obj. = %.2f + %.2f = %.2f\n',value(objP1),...
                value(objP2),value(objRob)));
        else
            title(sprintf('Obj. = %.2f - %.2f = %.2f\n',value(objP1),...
                abs(value(objP2)),value(objRob)));
        end

        saveas(gcf,sprintf('../Figs/FigRect_%s_%d.pdf',type,deg))
    end
    data.objPrim = objPrim;
    data.objPrimAg1 = objPrimAg1;
    data.objPrimAg2 = objPrimAg2;
    data.objDual = objDual;
end