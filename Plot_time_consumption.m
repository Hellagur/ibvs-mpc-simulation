function Plot_time_consumption(param, hist, k, saveFig)
%Plot_TIME_CONSUMPTION
%   Plots the computational time consumption of the MPC solver at each timestep.
%
% INPUTS:
%   param   - struct containing system parameters (param.ts: simulation timestep)
%   hist    - struct containing historical data, including:
%               hist.tcon - time consumption at each timestep
%   k       - total number of timesteps
%   saveFig - boolean flag to export the figure as PDF
%
% OUTPUT:
%   None (plots figure to MATLAB and optionally saves as PDF)

    % === Create figure ===
    fig = figure('Units','inches','Position',[1 1 8 6]);

    % === Time vector for x-axis ===
    tspan = (0:1:k) * param.ts;  % time vector in seconds

    % === Plot time consumption of each timestep ===
    % hist.tcon should have length k
    plot(tspan(1:end-1), hist.tcon, 'b-', 'LineWidth', 1.5);  

    % === Axis labels and formatting ===
    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Time Consumption (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
    grid on;                  % show grid
    axis tight;               % fit axes tightly to data
    box on;                   % draw box around plot
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');  % set font properties

    %% === Export figure as PDF if requested ===
    if saveFig == true
        set(gcf, 'PaperPositionMode', 'auto');  % maintain figure size
        set(gcf, 'Renderer', 'painters');       % use vector graphics for PDF
        figure_name = strcat('time_consumption_duration=', num2str(k), 's');
        print(fig, figure_name, '-dpdf', '-r600');  % save high-resolution PDF
    end
end
