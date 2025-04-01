function plotSystemResponse(s11, t, plotPos, filename)
    % Create figure and set position
    fig = figure();
    fig.Position = plotPos;
    
    % Create tiled layout
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Plot system output
    nexttile;
    plot(t, s11.y, 'b', 'LineWidth', 1.5);
    title('System Output y_t');
    xlabel('Time Step');
    ylabel('y_t');
    ylim padded;
    xlim padded;
    grid on;
    
    % Plot control input
    nexttile;
    plot(t, s11.u, 'r', 'LineWidth', 1.5);
    title('Control Input u_t');
    xlabel('Time Step');
    ylabel('u_t');
    ylim padded;
    xlim padded;
    grid on;
    
    % Plot reference omega
    nexttile;
    plot(t, s11.omega, 'b', 'LineWidth', 1.5);
    title('Ref \omega');
    xlabel('Time Step');
    ylabel('\omega');
    ylim padded;
    xlim padded;
    grid on;
    
    % Plot error
    nexttile;
    plot(t, s11.omega - s11.y, 'r', 'LineWidth', 1.5);
    title('Error');
    xlabel('Time Step');
    ylabel('error');
    ylim padded;
    xlim padded;
    grid on;
    
    % Save the figure
    filePath = fullfile("output", filename);
    exportgraphics(fig, filePath, 'Resolution', 300);
end