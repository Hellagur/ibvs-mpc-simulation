function Plot_relative_attitude_trajectory(param, hist, k, saveFig)
%PLOT_RELATIVE_ATTITUDE_TRAJECTORY
%   Plot the relative attitude trajectory (MRPs) of the chaser w.r.t target.
%   - Computes relative attitude in MRPs from DCMs
%   - Plots σ1, σ2, σ3 over time
%   - Adds zoomed inset for the last 50 seconds
%   - Optionally saves figure as PDF
%
% Inputs:
%   param    - struct with fields:
%                ts : sampling time
%   hist     - struct with fields:
%                R_ci : cell array of chaser DCMs (3x3)
%                R_ti : cell array of target DCMs (3x3)
%   k        - number of steps
%   saveFig  - logical, whether to save figure

    % Create figure
    fig = figure('Units','inches','Position',[1 1 8 3]);

    % Time vector
    tspan = (0:1:k)*param.ts;

    % Preallocate MRPs of relative attitude (σ_CT)
    s_ct  = zeros(3, k+1);

    %% ---- Compute relative attitude MRPs ----
    for i = 1:k+1
        R_ct = hist.R_ci{i} * hist.R_ti{i}';  % Relative DCM: chaser w.r.t target
        s_ct(:,i) = dcm2mrp(R_ct);            % Convert DCM to MRP
    end

    %% ---- Determine zoom-in range (last 50 seconds) ----
    zoom_duration = 50.0;  % seconds
    [~, zoom_start_idx] = min(abs(tspan - (tspan(end) - zoom_duration)));
    zoom_t = tspan(zoom_start_idx:end);

    %% ---- Plot MRPs ----
    hold on; box on;
    plot(tspan, s_ct(1,:), 'r-', 'LineWidth', 1.5);  % σ1
    plot(tspan, s_ct(2,:), 'k-', 'LineWidth', 1.5);  % σ2
    plot(tspan, s_ct(3,:), 'b-', 'LineWidth', 1.5);  % σ3

    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Relative Attitude (MRPs)', 'FontSize', 12, 'FontName', 'Times New Roman');
    grid on; axis tight;
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

    legend({'$\sigma_{CT,1}$', '$\sigma_{CT,2}$', '$\sigma_{CT,3}$'}, ...
           'Interpreter','latex', 'FontSize', 10, 'Location','northeast');

    %% ---- Inset: zoomed MRPs ----
    ax_inset = axes('Position', [0.65, 0.40, 0.20, 0.20]); hold on; box on;
    plot(zoom_t, s_ct(1,zoom_start_idx:end), 'r-', 'LineWidth', 1.5);
    plot(zoom_t, s_ct(2,zoom_start_idx:end), 'k-', 'LineWidth', 1.5);
    plot(zoom_t, s_ct(3,zoom_start_idx:end), 'b-', 'LineWidth', 1.5);
    set(gca,'FontSize',10,'FontName','Times New Roman'); 
    grid on; axis tight;
    ax_inset.YAxis.Exponent = 0;  % avoid scientific notation

    %% ---- Export figure as PDF if requested ----
    if saveFig
        set(gcf, 'PaperPositionMode', 'auto'); % Auto paper size
        set(gcf, 'Renderer', 'painters');      % Vector renderer
        figure_name = strcat('figs/relative_attitude_trajectory_duration=', ...
                             num2str(k*param.ts), 's');
        print(fig, figure_name, '-dpdf', '-r600');
    end
end
