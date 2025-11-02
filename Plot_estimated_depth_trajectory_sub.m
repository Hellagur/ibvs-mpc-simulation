function Plot_estimated_depth_trajectory_sub(~, hist, k, saveFig)
%PLOT_ESTIMATED_DEPTH_TRAJECTORY_SUB
%   Plots the estimated depth trajectories of feature points over time.
%   Each subplot corresponds to one feature point.
%
% Inputs:
%   ~        — unused (placeholder for consistency)
%   hist     — Struct containing:
%                 hist.rc   : 3N × k true feature coordinates in camera frame
%                 hist.zhat : N × k estimated depths
%   k        — Scalar, current time step index
%   saveFig  — Logical, true to save figure as PDF

    fig = figure('Units','inches','Position',[1 1 8 6]);

    N = size(hist.zhat, 1);  % number of feature points
    for i = 1:N
        subplot(2,2,i);
        % Plot true Z-coordinate (red solid line)
        plot(hist.rc(3*i,1:k),'r', 'LineWidth', 1.5); hold on;
        % Plot estimated depth (blue dashed line)
        plot(hist.zhat(i,1:k),'b--', 'LineWidth', 1.5);
        grid on; axis tight;

        xlabel('Time (s)', 'FontSize', 10, 'FontName', 'Times New Roman');
        ylabel(['$$\rm{Depth}\ Z_{' num2str(i) '}$$'], 'FontSize', 10, ...
            'FontName', 'Times New Roman', 'Interpreter', 'latex');

        ylim([0,22]); yticks(0:2:22);
        set(gca,'FontSize',10,'FontName','Times New Roman');
    end

    %% Save figure if requested
    if saveFig
        set(gcf, 'PaperPositionMode', 'auto');  % automatic paper size
        set(gcf, 'Renderer', 'painters');       % vector graphics
        % param.version needs to be in workspace or passed as input
        figure_name = strcat('plot_estimated_depth_k=', num2str(k));  
        print(fig, figure_name, '-dpdf', '-r600');
    end
end
