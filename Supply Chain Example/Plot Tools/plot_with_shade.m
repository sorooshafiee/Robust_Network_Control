function p = plot_with_shade(x, Y, plt)

    fig = figure;
    set(fig, 'Units', 'normalized', 'Position', plt.position)
    
    num = size(Y, 3);
    p = cell(num, 1);
    hold on
    box on
    for i = 1 : num
        p{i} = plot(x, mean(Y(:,:,i), 2), 'linewidth', 3, 'color', plt.colors(i, :));
        x2 = [x, flip(x)];
        fill(x2,[prctile(Y(:,:,i), plt.prc, 2)', flip(prctile(Y(:,:,i), 100-plt.prc, 2))'], ..., 
            plt.colors(i, :), 'LineStyle', 'none');
        alpha(plt.alpha)
    end
    
    set(gca, 'FontSize', plt.font_size - 6, 'TickLabelInterpreter', 'latex');
    xlabel(plt.xlabel, 'Interpreter', 'latex', 'FontSize', plt.font_size);
    ylabel(plt.ylabel, 'Interpreter', 'latex', 'FontSize', plt.font_size);
    %ylim([-5,120])
    
    if plt.grid
        grid on
    end
    
    if ~isempty(plt.lgd)
        lgd = legend([p{:}], plt.lgd, 'Interpreter', 'latex', 'Location', 'southeast');
        lgd.FontSize = plt.font_size - 3;
    end
    
    remove_border()
    saveas(gcf, plt.path, 'svg')
end