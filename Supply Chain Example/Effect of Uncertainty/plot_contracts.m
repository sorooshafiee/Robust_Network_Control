function plot_contracts(b_l, b_u, T, plt)
    fig = figure;
    set(fig, 'Units', 'normalized', 'Position', [0.35, 0.25, 0.4, 0.55])
    hold on
    grid on
    box on
    plot(1:T, b_u(1,:), '--', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    plot(1:T, b_u(2,:), 'color', [0.4660, 0.6740, 0.1880], 'linewidth', 2);
    plot(1:T, b_u(3,:), '.-.', 'color', [0.6350, 0.0780, 0.1840], 'linewidth', 2);
    
    plot(1:T, b_l(1,:), '--', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    plot(1:T, b_l(2,:), 'color', [0.4660, 0.6740, 0.1880], 'linewidth', 2);
    plot(1:T, b_l(3,:), '.-.', 'color', [0.6350, 0.0780, 0.1840], 'linewidth', 2);
    
    set(gca, 'FontSize', plt.font_size - 6, 'TickLabelInterpreter', 'latex');
    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', plt.font_size);
    ylabel('QF bounds', 'Interpreter', 'latex', 'FontSize', plt.font_size)
    title(plt.title, 'Interpreter', 'latex', 'FontSize', plt.font_size)
    lgd = legend('$\theta = 0.25$', '$\theta = 0.5$' , '$\theta = 1$', ...
                 'Interpreter', 'latex', 'Location', 'southwest');
    lgd.FontSize = plt.font_size - 3;
    ylim([0.5,5.5]);
    xlim([0, T+1]);
    remove_border()
    saveas(gcf, plt.path, 'svg')
end