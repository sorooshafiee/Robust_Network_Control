function plot_objective_functions(dataCent, dataDecent, data)

dataMat = zeros(37,7);
for kk = 1:37
    dataMat(kk,1) = dataCent.objPrim;
    dataMat(kk,2) = dataDecent.objPrim;
    dataMat(kk,3) = data{kk,1}.objPrim;
    dataMat(kk,4) = data{kk,2}.objPrim;
    dataMat(kk,5) = data{kk,3}.objPrim;
    dataMat(kk,6) = data{kk,4}.objPrim;
end


% create the figures
tV = 0:5:180;



% Set up the figure properties
fig = figure;
fig.Color = 'w';
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.Units = 'inches';
fig.PaperPosition = [0, 0, 5, 3.09];
fig.PaperSize = [5, 3.09];
fig.Position = [0.1, 0.1, 4.9, 2.99];
fig.Resize = 'off';
fig.InvertHardcopy = 'off';


h = zeros(6,1);
ax = gca;
hold on;
grid on
idx = find(dataMat(:,1)<=40);
h(1) = plot(tV(idx),dataMat(idx,1), 'color', [0, 0.4470, 0.7410]);
h(2) = plot(tV(idx),dataMat(idx,2),'k--');
h(3) = plot(tV(idx),dataMat(idx,3),'-s', 'color', [0.4660, 0.6740, 0.1880], 'MarkerSize', 4);
h(4) = plot(tV,dataMat(:,4),'-o', 'color', [0.6350, 0.0780, 0.1840], 'MarkerSize', 4);

idx = find(dataMat(:,5)>=40);
h(5) = plot(tV(1:idx(1)-1),dataMat(1:idx(1)-1,5),'-d', 'color', [0.6350, 0.0780, 0.1840], 'MarkerSize', 4);
plot(tV(idx(end)+1:end),dataMat(idx(end)+1:end,5),'-d', 'color', [0.6350, 0.0780, 0.1840], 'MarkerSize', 4);

idx = find(dataMat(:,6)>=40);
h(6) = plot(tV(1:idx(1)-1),dataMat(1:idx(1)-1,6),'-*', 'color', [0.6350, 0.0780, 0.1840], 'MarkerSize', 4);
plot(tV(idx(end)+1:end),dataMat(idx(end)+1:end,6),'-*', 'color', [0.6350, 0.0780, 0.1840], 'MarkerSize', 4);



ax.FontName = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.Title.Interpreter = 'LaTeX';
ax.XLabel.Interpreter = 'LaTeX';
ax.YLabel.Interpreter = 'LaTeX';
ax.Box = 'on';
% ax.LineWidth = 1.5;
ax.FontSize = 8;

ax.XLim = [0, 180];
ax.YLim = [-3, 30];
ax.XTick = 0:20:180;    % The locations of the tick marks
ax.YTick = 0:5:30;



xlabel('Degrees');
ylabel('Objective value');


legend(h,{'$\mathcal{CP}$','$\mathcal{NS}$','$\mathcal{DP}$ - Rect.',...
    '$\mathcal{DP}$ - El. sc.=1.5','$\mathcal{DP}$ - El. sc.=3','$\mathcal{DP}$ - El. sc.=10'},...
    'Position',[0.25 0.7 0.1 0.1],'Interpreter','Latex');


saveas(gcf,'../Figs/CostPlot.pdf')