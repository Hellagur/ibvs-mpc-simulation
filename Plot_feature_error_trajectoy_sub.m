function Plot_feature_error_trajectoy_sub(param, hist, k, saveFig)
%PLOT_FEATURE_ERROR_TRAJECTOY_SUB
%   Plot subplots of each feature error component (u/v) over time.
%   Each subplot includes a main curve and an inset zoom for the last 50 seconds.
%
% Inputs:
%   param    — Struct containing simulation parameters (ts, sd, version)
%   hist     — Struct containing simulation history (xs)
%   k        — Current simulation step
%   saveFig  — Logical, true to save figure as PDF
%
% Features:
%   - 4 features, 2 components each (u,v) → 8 subplots
%   - Main plot + inset zoom for terminal error
%   - LaTeX labels, Times New Roman font

    fig = figure('Units','inches','Position',[1 1 12 9]);
    tspan = (0:1:k) * param.ts;  % simulation time vector
    
    for i = 1:8
        %% ---- Create main subplot ----
        ax_main = subplot(4,2,i); 
        view(2); hold(ax_main, 'on');
        
        % Determine error component (u or v) and feature index
        sub_prefix = {'u', 'v'};
        err = hist.xs(i,1:k+1) - param.sd(i);  % deviation from desired
        
        % Plot main error trajectory
        plot(tspan, err, 'bo', 'LineWidth', 0.5);
        grid on; axis tight;
        
        % Label with LaTeX
        sub_name = strcat('$', sub_prefix{mod(i-1,2)+1}, '_', num2str(floor((i-1)/2)+1), '\ (\rm{pixel})$');
        xlabel('Time (s)', 'FontSize', 10, 'FontName', 'Times New Roman');
        ylabel(sub_name, 'FontSize', 10, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
        set(ax_main,'FontSize',10,'FontName','Times New Roman');

        %% ---- Inset for terminal error (last 50s) ----
        zoom_duration = 50.0;  % seconds
        [~, zoom_start_idx] = min(abs(tspan - (tspan(end) - zoom_duration)));
        zoom_t = tspan(zoom_start_idx:end);
        zoom_err = err(zoom_start_idx:end);
    
        % Get main subplot position
        pos = get(ax_main, 'Position');  % [x y width height]
    
        % Set inset size (relative to main plot)
        inset_scale = 0.27;  
        inset_w = pos(3) * inset_scale;
        inset_h = pos(4) * inset_scale;
    
        % Position inset in upper right corner of main subplot
        inset_x = pos(1) + pos(3) - inset_w - 0.01;
        inset_y = pos(2) + pos(4) - inset_h - 0.01;
    
        % Create inset axes
        ax_inset = axes('Position', [inset_x, inset_y, inset_w, inset_h]);
        box on;
        plot(ax_inset, zoom_t, zoom_err, 'r', 'LineWidth', 1);
        title(ax_inset, '', 'FontSize', 8);
        set(ax_inset, 'FontSize', 8, 'FontName', 'Times New Roman');
        grid(ax_inset, 'on');
        axis(ax_inset, 'tight');
        % xticks can be customized if needed
        % xticks(ax_inset, 96:2:100)
    end
    
    %% ---- Save figure as PDF ----
    if saveFig
        set(gcf, 'PaperPositionMode', 'auto');  % Auto paper size
        set(gcf, 'Renderer', 'painters');       % Vector rendering
        figure_name = strcat('plot_feature_trajectory2D_k=', num2str(k), param.version);
        print(fig, figure_name, '-dpdf', '-r600');
    end
end
