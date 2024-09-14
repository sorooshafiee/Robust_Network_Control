function data = toyExampleDecent(con1,con2,objP1,objP2,conW1,conW2,conFS,x1,x2,u1,u2,w1,w2,fDual,fPlot)

    % constraints
    conP = con1 + con2;
    objP = objP1 + objP2;

    % add constaint of decentalized problem
    conP = conP + [u1 == sdpvar(1,1) + sdpvar(1,2)*w1];

    % uncertainty
    conW = conW1 + conW2;

    % formulate robust optimization problem
    x = [x1;x2];
    u = [u1;u2];
    w = [w1;w2];
    wE = [];
    uLD = [];


    robustP = RobustProblem(objP,conP,conW,x,u,uLD,[],w,wE);
    robustP.robustify(1);
    robustP.getRobust;

    objRob = robustP.robustPB.object;
    conRob = robustP.robustPB.constr;

    optimize(conRob,objRob,sdpsettings('verbose',0));
    fprintf('Decentralized Problem: Robust problem objective: %f + %f = %f\n',value(objP1),...
        value(objP2),value(objRob));

    objPrim = value(objRob);
    objPrimAg1 = value(objP1);
    objPrimAg2 = value(objP2);

    if fDual
        robustP.scenario;
        robustP.getScenarioLB;
        fprintf('Centralized Problem: ');
        objDual = robustP.solveScenarioLB;
    else
        objDual = NaN;
    end


    if fPlot
        % plot the feasible region

        x110 = value(robustP.robustPB.x{1,2}(1));
        x111 = value(robustP.robustPB.x{1,2}(2));
        x112 = value(robustP.robustPB.x{1,2}(3));
        x113 = value(robustP.robustPB.x{1,2}(4));
        x114 = value(robustP.robustPB.x{1,2}(5));


        x120 = value(robustP.robustPB.x{2,2}(1));
        x121 = value(robustP.robustPB.x{2,2}(2));
        x122 = value(robustP.robustPB.x{2,2}(3));
        x123 = value(robustP.robustPB.x{2,2}(4));
        x124 = value(robustP.robustPB.x{2,2}(5));



        conf = conW1 + conW2 + [x1(1,2) == x110 + [x111,x112]*w1 + [x113,x114]*w2];
        conf = conf + [x1(2,2) == x120 + [x121,x122]*w1 + [x123,x124]*w2];

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
        plot(conFS,x1(:,2),[244, 182, 182]/255);
        plot(conf,x1(:,2),[179, 179, 249]/255);
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
        saveas(gcf,'../Figs/FigDecent.pdf')
    end
    data.objPrim = objPrim;
    data.objPrimAg1 = objPrimAg1;
    data.objPrimAg2 = objPrimAg2;
    data.objDual = objDual;
end