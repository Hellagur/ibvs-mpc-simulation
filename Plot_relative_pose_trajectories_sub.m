function Plot_relative_pose_trajectories_sub(param, hist, k, saveFig)
%PLOT_RELATIVE_POSE_TRAJECTORIES_SUB
%   Plot relative position and attitude (MRPs) trajectories in subplots.
%   Includes zoom-in insets for the last 50 seconds.
%
% Inputs:
%   param    - struct with sampling time and limits
%   hist     - struct with fields:
%                R_ti, R_li, R_ci : rotation matrices (cell arrays)
%                xl              : relative position history
%   k        - number of steps
%   saveFig  - logical, whether to save the figure as PDF

    % Create figure
    fig = figure('Units','inches','Position',[1 1 8 6]);

    % Time vector
    tspan = (0:1:k) * param.ts;

    % Preallocate relative position and attitude arrays
    r_ct_t = zeros(3, k+1);  % Relative position in target frame
    s_ct   = zeros(3, k+1);  % Relative attitude (MRPs)

    %% ---- Compute relative position ----
    for i = 1:k+1
        R_tl = hist.R_ti{i} * hist.R_li{i}';     % DCM from LVLH to target frame
        r_ct_t(:,i) = R_tl * hist.xl(7:9,i);
    end

    %% ---- Compute relative attitude (MRPs) ----
    for i = 1:k+1
        R_ct = hist.R_ci{i} * hist.R_ti{i}';     % DCM from target to chaser
        s_ct(:,i) = dcm2mrp(R_ct);               % Convert DCM to MRP
    end

    %% ---- Determine zoom-in range (last 50 seconds) ----
    zoom_duration = 50.0;
    [~, zoom_start_idx] = min(abs(tspan - (tspan(end) - zoom_duration)));
    zoom_t = tspan(zoom_start_idx:end);

    %% ---- Subplot 1: Relative position ----
    subplot(211);
    hold on; box on;
    plot(tspan, r_ct_t(1,:), 'r-', 'LineWidth', 1.5);
    plot(tspan, r_ct_t(2,:), 'k-', 'LineWidth', 1.5);
    plot(tspan, r_ct_t(3,:), 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Rel. Position in $\mathcal{F}_T$ (m)', 'FontSize', 12, ...
           'FontName', 'Times New Roman', 'Interpreter','latex');
    grid on; axis tight;
    ylim([min(r_ct_t,[],'all')-1, max(r_ct_t,[],'all')+1]);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
    legend({'$x_T$', '$y_T$', '$z_T$'}, 'Interpreter','latex', ...
           'FontSize', 10, 'Location','northeast');

    % Inset for x,y
    ax_inset_xy = axes('Position', [0.38, 0.63, 0.20, 0.10]); hold on; box on;
    plot(zoom_t, r_ct_t(1,zoom_start_idx:end), 'r-', 'LineWidth', 1.5);
    plot(zoom_t, r_ct_t(2,zoom_start_idx:end), 'k-', 'LineWidth', 1.5);
    set(gca,'FontSize',10,'FontName','Times New Roman');
    ax_inset_xy.YAxis.Exponent = 0;  % Disable scientific notation
    grid on; axis tight;

    % Inset for z
    ax_inset_z = axes('Position', [0.65, 0.63, 0.20, 0.10]); hold on; box on;
    plot(zoom_t, r_ct_t(3,zoom_start_idx:end), 'b-', 'LineWidth', 1.5);
    set(gca,'FontSize',10,'FontName','Times New Roman');
    ax_inset_z.YAxis.Exponent = 0;
    grid on; axis tight;

    %% ---- Subplot 2: Relative attitude (MRPs) ----
    subplot(212);
    hold on; box on;
    plot(tspan, s_ct(1,:), 'r-', 'LineWidth', 1.5);
    plot(tspan, s_ct(2,:), 'k-', 'LineWidth', 1.5);
    plot(tspan, s_ct(3,:), 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Rel. Attitude $\mathbf{\sigma}_{CT}$ (MRPs)', 'FontSize', 12, ...
           'FontName', 'Times New Roman', 'Interpreter','latex');
    grid on; axis tight;
    ylim([min(s_ct,[],'all')-0.01, max(s_ct,[],'all')+0.01]);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
    legend({'$\sigma_{CT,1}$', '$\sigma_{CT,2}$', '$\sigma_{CT,3}$'}, ...
           'Interpreter','latex', 'FontSize', 10, 'Location','northeast');

    % Inset for attitude (σ1, σ2, σ3)
    ax_inset_att = axes('Position', [0.65, 0.22, 0.20, 0.10]); hold on; box on;
    plot(zoom_t, s_ct(1,zoom_start_idx:end), 'r-', 'LineWidth', 1.5);
    plot(zoom_t, s_ct(2,zoom_start_idx:end), 'k-', 'LineWidth', 1.5);
    plot(zoom_t, s_ct(3,zoom_start_idx:end), 'b-', 'LineWidth', 1.5);
    set(gca,'FontSize',10,'FontName','Times New Roman');
    ax_inset_att.YAxis.Exponent = 0;
    grid on; axis tight;

    %% ---- Export figure as PDF if requested ----
    if saveFig
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Renderer', 'painters');
        figure_name = strcat('figs/relative_pose_trajectories_duration=', ...
                             num2str(k*param.ts), 's');
        print(fig, figure_name, '-dpdf', '-r600');
    end
end
