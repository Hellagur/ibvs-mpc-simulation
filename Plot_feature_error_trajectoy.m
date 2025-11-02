function Plot_feature_error_trajectoy(param, hist, k, saveFig)
%PLOT_FEATURE_ERROR_TRAJECTOY
%   Plot the 2D feature position errors (u, v) in pixels over time.
%   Includes gradient coloring along time and an inset zoom for the terminal period.
%
% Inputs:
%   param    — Struct containing simulation parameters (e.g., ts, sd)
%   hist     — Struct containing simulation history (e.g., xs for feature positions)
%   k        — Current simulation time step
%   saveFig  — Logical, true to save figure as PDF
%
% Features:
%   - Gradient coloring for each feature over time
%   - Separate colors for u and v components
%   - Main plot + inset for last 50s
%   - Custom legend and axis formatting

    fig = figure('Units','inches','Position',[1 1 8 6]);
    
    tspan = (0:1:k)*param.ts;
    view(2); hold on;
    
    %% ---- Base colors for features (u and v components) ----
    base_colors = [... % RGB values
        231, 76,  60;    % feature1 u
        213, 32,  60;    % feature1 v
        46,  204, 113;   % feature2 u
        46,  172, 178;   % feature2 v
        52,  152, 219;   % feature3 u
        132, 136, 235;   % feature3 v
        241, 196, 15;    % feature4 u
        241, 132, 15];   % feature4 v
    base_colors = base_colors / 255;
    
    markerSize = 28;  % scatter marker size
    
    %% ---- Compute error ----
    err = hist.xs(1:8,1:k+1) - param.sd;  % deviation from desired
    
    h_legend = gobjects(8,1);  % legend handles
        
    %% ---- Plot each feature trajectory ----
    for feat = 1:4
        idx_u = 2*feat-1; % u index
        idx_v = 2*feat;   % v index
        
        color_u_base = base_colors(idx_u,:);
        color_v_base = base_colors(idx_v,:);
        
        % u trajectory
        for j = 1:length(tspan)
            alpha = 0.5 + 0.5*(j-1)/k;  % gradient intensity
            color = color_u_base * alpha + (1-alpha)*[1 1 1];
            scatter(tspan(j), err(idx_u,j), markerSize, color, 'filled', ...
                'MarkerEdgeColor', color_u_base, 'LineWidth', 0.75, ...
                'MarkerFaceAlpha',0.5); hold on;
            if j > 1
                plot(tspan(j-1:j), err(idx_u,j-1:j), '-', 'Color', color, 'LineWidth', 1.2);
            end
        end
        
        % v trajectory (same as u)
        for j = 1:length(tspan)
            alpha = 0.5 + 0.5*(j-1)/k;
            color = color_v_base * alpha + (1-alpha)*[1 1 1];
            scatter(tspan(j), err(idx_v,j), markerSize, color, 'filled', ...
                'MarkerEdgeColor', color_v_base, 'LineWidth', 0.75, ...
                'MarkerFaceAlpha',0.5); hold on;
            if j > 1
                plot(tspan(j-1:j), err(idx_v,j-1:j), '-', 'Color', color, 'LineWidth', 1.2);
            end
        end
        
        %% ---- Legend entries ----
        h_u(feat) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', color_u_base, ...
            'MarkerFaceColor', color_u_base, 'MarkerSize', 6, 'DisplayName', sprintf('u_%d', feat));
        h_v(feat) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', color_v_base, ...
            'MarkerFaceColor', color_v_base, 'MarkerSize', 6, 'DisplayName', sprintf('v_%d', feat));
        h_legend(feat)   = h_u(feat);
        h_legend(feat+4) = h_v(feat);
    end
    
    %% ---- Axis formatting ----
    grid on; axis tight;
    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Position Error (pixel)', 'FontSize', 12, 'FontName', 'Times New Roman');
    set(gca,'FontSize',12,'FontName','Times New Roman');
    
    %% ---- Legend ----
    lgd = legend(h_legend, ...
        {'$u_1$','$u_2$','$u_3$','$u_4$', ...
         '$v_1$','$v_2$','$v_3$','$v_4$'}, ...
        'NumColumns',2, 'FontSize',10, 'FontName','Times New Roman', ...
        'Location','northeast','Interpreter','latex');
    lgd.Box = 'on';
    lgd.ItemTokenSize = [18,9];
    
    %% ---- Inset for terminal error zoom (last 50s) ----
    zoom_duration = 50.0;  % seconds
    [~, zoom_start_idx] = min(abs(tspan - (tspan(end) - zoom_duration)));
    zoom_t = tspan(zoom_start_idx:end);
    zoom_err = err(:,zoom_start_idx:end);
    
    pos = get(gca, 'Position');
    inset_scale = 0.35;
    ax_inset = axes('Position', [pos(1)+pos(3)-pos(3)*inset_scale-0.05, ...
                                 pos(2)+pos(4)-pos(4)*inset_scale-0.45, ...
                                 pos(3)*inset_scale, pos(4)*inset_scale]);
    box on; hold on;
    
    for feat = 1:4
        idx_u = 2*feat-1; idx_v = 2*feat;
        color_u_base = base_colors(idx_u,:);
        color_v_base = base_colors(idx_v,:);
        for j = 1:length(zoom_t)
            alpha = 0.5 + 0.5*(j-1)/(length(zoom_t)-1);
            color_u = color_u_base * alpha + (1-alpha)*[1 1 1];
            color_v = color_v_base * alpha + (1-alpha)*[1 1 1];
            scatter(ax_inset, zoom_t(j), zoom_err(idx_u,j), markerSize, color_u, 'filled', ...
                'MarkerEdgeColor', color_u_base, 'LineWidth', 0.75, 'MarkerFaceAlpha',0.5);
            scatter(ax_inset, zoom_t(j), zoom_err(idx_v,j), markerSize, color_v, 'filled', ...
                'MarkerEdgeColor', color_v_base, 'LineWidth', 0.75, 'MarkerFaceAlpha',0.5);
            if j > 1
                plot(ax_inset, zoom_t(j-1:j), zoom_err(idx_u,j-1:j), '-', 'Color', color_u, 'LineWidth', 1.2);
                plot(ax_inset, zoom_t(j-1:j), zoom_err(idx_v,j-1:j), '-', 'Color', color_v, 'LineWidth', 1.2);
            end
        end
    end
    set(ax_inset, 'FontSize',10,'FontName','Times New Roman');
    grid(ax_inset,'on'); axis(ax_inset,'tight');
    
    %% ---- Save figure as PDF ----
    if saveFig
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Renderer', 'painters');
        figure_name = strcat('figs/feature_error_trajectory_duration=', num2str(k*param.ts), 's');
        print(fig, figure_name, '-dpdf', '-r600');
    end
end
