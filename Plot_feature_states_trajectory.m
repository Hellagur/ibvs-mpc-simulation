function Plot_feature_states_trajectory(param, hist, k, saveFig)
%PLOT_FEATURE_STATES_TRAJECTORY
%   Plot the trajectories of feature points in the image plane (u,v).
%   - Draws background grid and axes
%   - Plots trajectories of 4 feature points with color gradient to show time progression
%   - Marks initial, final, target positions
%   - Optionally saves the figure as PDF
%
% Inputs:
%   param    - struct with parameters:
%                um, nm : image plane range
%                sd     : desired feature positions
%                ts     : sampling period
%   hist     - struct with history data:
%                xs     : feature states (8 x k+1)
%   k        - current step index
%   saveFig  - logical, whether to save figure

    fig = figure('Units','inches','Position',[1 1 8 6]);

    %% ---- Draw background grid ----
    step = 20;
    x_range = -param.um:step:param.um;
    y_range = -param.nm:step:param.nm;
    hold on;
    
    % Horizontal grid lines
    for y = y_range
        line([-param.um, param.um], [y, y], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.6);
    end
    % Vertical grid lines
    for x = x_range
        line([x, x], [-param.nm, param.nm], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.6);
    end
    
    % Background color and box
    set(gca, 'Color', [0.99 0.99 0.99]);  
    box on;
    
    % Center lines
    line([0 0], [-param.nm, param.nm], 'Color', [0.4 0.4 0.4], 'LineStyle', ':');
    line([-param.um, param.um], [0 0], 'Color', [0.4 0.4 0.4], 'LineStyle', ':');

    % Desired positions connected by dashed lines
    line('XData', [param.sd(1:2:7);param.sd(1)], ...
         'YData', [param.sd(2:2:8);param.sd(2)], ...
         'LineStyle','--','LineWidth',1.5);

    %% ---- Plot feature trajectories with color gradient ----
    base_colors = [...
        231, 76,  60;    % Red
        46,  204, 113;   % Green
        52,  152, 219;   % Blue
        241, 196, 15];   % Yellow
    base_colors = base_colors / 255;  % Normalize to [0,1]

    markerSize = 50; % Marker size
    for idx = 1:4
        xs_idx = hist.xs(2*idx-1,:); % X positions
        ys_idx = hist.xs(2*idx,:);   % Y positions
        % Linear gradient: light color → base color
        for i = 1:k
            color = 0.3 + 0.7*(i-1)/(k-1);   % Start at 0.3 (light)
            c_rgb = base_colors(idx,:) * color + (1-color)*[1 1 1]; 
            scatter(xs_idx(i), ys_idx(i), markerSize, c_rgb, 'filled'); hold on;
        end
        % Plot trajectory line
        plot(xs_idx, ys_idx, '-', 'Color', base_colors(idx,:), 'LineWidth', 1.0);
    end

    %% ---- Plot desired positions ----
    for i = 1:4
        p0 = plot(param.sd(2*i-1), param.sd(2*i), 'p', ...
             'MarkerFaceColor','#FF9966','MarkerEdgeColor','k','MarkerSize',11);
    end

    %% ---- Plot position at t = 70s (optional) ----
    for i = 1:4
        p1 = plot(hist.xs(2*i-1,141), hist.xs(2*i,141), 'o', 'MarkerFaceColor','none', ...
              'MarkerEdgeColor','b', 'MarkerSize',6, 'LineWidth',1.8);
    end

    %% ---- Plot initial positions ----
    for i = 1:4
        p2 = plot(hist.xs(2*i-1,1), hist.xs(2*i,1), 'o', ...
             'MarkerFaceColor','none','MarkerEdgeColor','r','MarkerSize',6,'LineWidth',1.5);
    end

    %% ---- Plot final positions ----
    for i = 1:4
        p3 = plot(hist.xs(2*i-1,k), hist.xs(2*i,k), 'o', ...
             'MarkerFaceColor','none','MarkerEdgeColor','g','MarkerSize',6,'LineWidth',1.5);
    end

    %% ---- Set axes ----
    axis equal;
    xlim([-param.um, param.um]);
    ylim([-param.nm, param.nm]);
    xlabel('$u\ (\rm{px})$', 'FontSize',12,'FontName','Times New Roman','Interpreter','latex');
    ylabel('$v\ (\rm{px})$', 'FontSize',12,'FontName','Times New Roman','Interpreter','latex');
    set(gca,'FontSize',12,'FontName','Times New Roman');

    %% ---- Legend ----
    legend([p3, p2, p1, p0], {'Final Position','Initial Position','t = 70s Position','Target Position'}, ...
           'FontSize',12,'FontName','Times New Roman','Location','northeast');

    %% ---- Save figure as PDF if requested ----
    if saveFig
        set(gcf, 'PaperPositionMode','auto');
        set(gcf, 'Renderer','painters');
        figure_name = strcat('figs/feature_states_trajectory_duration=', num2str(k*param.ts), 's');
        print(fig, figure_name, '-dpdf', '-r600');
    end
end
