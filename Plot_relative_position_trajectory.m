function Plot_relative_position_trajectory(param, hist, k, saveFig)
%PLOT_RELATIVE_POSITION_TRAJECTORY
%   Plot the relative position trajectory of the chaser in the target frame.
%   - Computes relative position using DCMs
%   - Plots x, y, z components over time
%   - Adds zoomed insets for the last 50 seconds
%   - Optionally saves the figure as a PDF
%
% Inputs:
%   param    - struct with fields:
%                ts : sampling time
%   hist     - struct with fields:
%                R_ti : cell array of target DCMs
%                R_li : cell array of LVLH DCMs
%                xl   : position vector history (at least 7:9 for relative)
%   k        - number of steps
%   saveFig  - logical, whether to save figure

    % Create figure
    fig = figure('Units','inches','Position',[1 1 8 3]);

    % Time vector
    tspan = (0:1:k) * param.ts;

    % Preallocate relative position array (in target frame)
    r_ct_t = zeros(3, k+1);

    %% ---- Compute relative position in target frame ----
    for i = 1:k+1
        R_tl = hist.R_ti{i} * hist.R_li{i}';       % DCM from LVLH to target
        r_ct_t(:,i) = R_tl * hist.xl(7:9,i);       % Relative position vector
    end

    %% ---- Determine zoom-in range (last 50 seconds) ----
    zoom_duration = 50.0;  % seconds
    [~, zoom_start_idx] = min(abs(tspan - (tspan(end) - zoom_duration)));
    zoom_t = tspan(zoom_start_idx:end);

    %% ---- Plot x, y, z over time ----
    hold on; box on;
    plot(tspan, r_ct_t(1,:), 'r-', 'LineWidth', 1.5);  % x_T
    plot(tspan, r_ct_t(2,:), 'k-', 'LineWidth', 1.5);  % y_T
    plot(tspan, r_ct_t(3,:), 'b-', 'LineWidth', 1.5);  % z_T

    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Rel. Pos in $\mathcal{F}_T$ (m)', 'FontSize', 12, ...
           'FontName', 'Times New Roman', 'Interpreter','latex');
    grid on; axis tight;
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

    legend({'$x_T$', '$y_T$', '$z_T$'}, 'Interpreter','latex', ...
        'FontSize', 10, 'Location','northeast');

    %% ---- Inset: zoom x and y ----
    ax_inset_xy = axes('Position', [0.38, 0.25, 0.20, 0.20]); hold on; box on;
    plot(zoom_t, r_ct_t(1,zoom_start_idx:end), 'r-', 'LineWidth', 1.5);
    plot(zoom_t, r_ct_t(2,zoom_start_idx:end), 'k-', 'LineWidth', 1.5);
    set(gca,'FontSize',10,'FontName','Times New Roman');
    ax_inset_xy.YAxis.Exponent = 0;  % Disable scientific notation
    grid on; axis tight;

    %% ---- Inset: zoom z ----
    ax_inset_z = axes('Position', [0.65, 0.25, 0.20, 0.20]); hold on; box on;
    plot(zoom_t, r_ct_t(3,zoom_start_idx:end), 'b-', 'LineWidth', 1.5);
    set(gca,'FontSize',10,'FontName','Times New Roman');
    ax_inset_z.YAxis.Exponent = 0;
    grid on; axis tight;

    %% ---- Export figure as PDF if requested ----
    if saveFig
        set(gcf, 'PaperPositionMode', 'auto');  % Auto paper size
        set(gcf, 'Renderer', 'painters');       % Vector renderer
        figure_name = strcat('figs/relative_position_trajectory_duration=', ...
                             num2str(k*param.ts), 's');
        print(fig, figure_name, '-dpdf', '-r600');
    end
end
