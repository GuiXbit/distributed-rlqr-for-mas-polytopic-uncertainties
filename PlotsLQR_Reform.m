function varargout = PlotsLQR_Reform(plotName, data, opts)
%PLOTSLQR_REFORM Centraliza os plots dos scripts nominal e robusto.
%
% Uso:
%   PlotsLQR_Reform(plotName, data)
%   PlotsLQR_Reform(plotName, data, opts)
%   h = PlotsLQR_Reform(plotName, data, opts)
%
% Entradas:
%   plotName : identificador do plot a ser gerado.
%   data     : struct com os dados necessarios para o plot.
%   opts     : struct opcional com configuracoes visuais.
%
% Esta funcao comeca como uma camada de padronizacao. Os plots serao
% adicionados aos poucos, conforme a definicao do layout desejado.

    if nargin < 2
        error('PlotsLQR_Reform requer ao menos plotName e data.');
    end

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    cfg = default_plot_options(opts);

    switch lower(string(plotName))
        case "state_trajectories"
            fig = plot_state_trajectories(data, cfg);
        case "epsilon_components"
            fig = plot_epsilon_components(data, cfg);
        case "epsilon_global_norm"
            fig = plot_epsilon_global_norm(data, cfg);
        case "control_inputs"
            fig = plot_control_inputs(data, cfg);
        case "phase_plane"
            fig = plot_phase_plane(data, cfg);
        case "phase_3d_time"
            fig = plot_phase_3d_time(data, cfg);
        case "uncertainty_norms"
            fig = plot_uncertainty_norms(data, cfg);
        case "costs"
            fig = plot_costs(data, cfg);
        case "tracking_global_and_components"
            fig = plot_tracking_global_and_components(data, cfg);
        otherwise
            error('Plot "%s" ainda nao foi implementado em PlotsLQR_Reform.', plotName);
    end

    if nargout > 0
        varargout{1} = fig;
    end
end

function cfg = default_plot_options(opts)
    cfg = struct();

    cfg.figure_name = get_option(opts, 'figure_name', '');
    cfg.figure_color = get_option(opts, 'figure_color', 'w');
    cfg.axes_font_size = get_option(opts, 'axes_font_size', 12);
    cfg.label_font_size = get_option(opts, 'label_font_size', 14);
    cfg.title_font_size = get_option(opts, 'title_font_size', 15);
    cfg.legend_font_size = get_option(opts, 'legend_font_size', 11);
    cfg.line_width_open = get_option(opts, 'line_width_open', 1.2);
    cfg.line_width_closed = get_option(opts, 'line_width_closed', 1.6);
    cfg.line_width_reference = get_option(opts, 'line_width_reference', 2.0);
    cfg.grid = get_option(opts, 'grid', 'on');
    cfg.box = get_option(opts, 'box', 'on');
    cfg.legend_location = get_option(opts, 'legend_location', 'best');
    cfg.interpreter = get_option(opts, 'interpreter', 'tex');
end

function fig = plot_state_trajectories(data, cfg)
    required_fields = {'t', 'x_open', 'x_closed', 'nstate', 'nagent'};
    validate_fields(data, required_fields, 'state_trajectories');
    has_leader = isfield(data, 'leader') && ~isempty(data.leader);

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);
    tl = tiledlayout(fig, data.nstate, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    legend_source_ax = gobjects(1,1);

    for idxState = 1:data.nstate
        ax_open = nexttile(tl);
        hold(ax_open, 'on');
        grid(ax_open, cfg.grid);
        box(ax_open, cfg.box);

        for i = 1:data.nagent
            plot(ax_open, data.t, squeeze(data.x_open(idxState,i,:)), '--', 'LineWidth', cfg.line_width_open);
        end
        if has_leader
            plot(ax_open, data.t, data.leader(idxState,:), 'k:', 'LineWidth', cfg.line_width_reference);
        end
        xlabel(ax_open, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        ylabel(ax_open, sprintf('x_%d', idxState), 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        title(ax_open, sprintf('x_%d^i(k) - Malha Aberta', idxState), 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
        ax_open.FontSize = cfg.axes_font_size;

        ax_closed = nexttile(tl);
        hold(ax_closed, 'on');
        grid(ax_closed, cfg.grid);
        box(ax_closed, cfg.box);

        for i = 1:data.nagent
            plot(ax_closed, data.t, squeeze(data.x_closed(idxState,i,:)), '-', 'LineWidth', cfg.line_width_closed);
        end
        if has_leader
            plot(ax_closed, data.t, data.leader(idxState,:), 'k:', 'LineWidth', cfg.line_width_reference);
        end
        xlabel(ax_closed, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        ylabel(ax_closed, sprintf('x_%d', idxState), 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        title(ax_closed, sprintf('x_%d^i(k) - Malha Fechada', idxState), 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
        ax_closed.FontSize = cfg.axes_font_size;

        if idxState == 1
            legend_source_ax = ax_open;
        end
    end

    legend_labels = get_state_legend_labels(data, has_leader);
    lgd = legend(flipud(legend_source_ax.Children), legend_labels, 'Location', 'southoutside', 'Orientation', 'horizontal');
    lgd.FontSize = cfg.legend_font_size;
    lgd.Interpreter = cfg.interpreter;
end

function fig = plot_epsilon_components(data, cfg)
    required_fields = {'t', 'epsilon', 'nstate', 'nagent'};
    validate_fields(data, required_fields, 'epsilon_components');

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);
    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, cfg.grid);
    box(ax, cfg.box);

    labels = cell(1, data.nstate * data.nagent);
    idxLabel = 1;
    for i = 1:data.nagent
        for j = 1:data.nstate
            plot(ax, data.t, squeeze(data.epsilon(j,i,:)), '-', 'LineWidth', cfg.line_width_closed);
            labels{idxLabel} = sprintf('\\epsilon_{%d%d}', i, j);
            idxLabel = idxLabel + 1;
        end
    end

    xlabel(ax, 'Iteration Steps', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ylabel(ax, 'Tracking Error Dynamics', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    title(ax, 'Tracking Error Dynamics', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
    ax.FontSize = cfg.axes_font_size;

    lgd = legend(ax, labels, 'Location', cfg.legend_location);
    lgd.FontSize = cfg.legend_font_size;
    lgd.Interpreter = cfg.interpreter;
end

function fig = plot_epsilon_global_norm(data, cfg)
    required_fields = {'t', 'epsilon'};
    validate_fields(data, required_fields, 'epsilon_global_norm');

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);
    nsteps = numel(data.t);
    eps_global_norm = zeros(1, nsteps);
    for k = 1:nsteps
        eps_global_norm(k) = norm(data.epsilon(:,:,k), 'fro');
    end

    has_tracking = isfield(data, 'x_closed') && isfield(data, 'leader') && ~isempty(data.leader);
    has_bound = has_tracking && isfield(data, 'LB') && ~isempty(data.LB);

    if has_tracking
        ax1 = subplot(1,2,1); hold(ax1, 'on');
        grid(ax1, cfg.grid);
        box(ax1, cfg.box);

        tracking_global = zeros(1, nsteps);
        for k = 1:nsteps
            tracking_mat = data.x_closed(:,:,k) - data.leader(:,k);
            tracking_global(k) = norm(tracking_mat(:), 2);
        end

        plot(ax1, data.t, tracking_global, '-', 'LineWidth', cfg.line_width_closed);
        if has_bound
            sigma_min_LB = min(svd(data.LB));
            tracking_bound = eps_global_norm / sigma_min_LB;
            plot(ax1, data.t, tracking_bound, 'k--', 'LineWidth', cfg.line_width_reference);
            lgd1 = legend(ax1, {'||[x^1-x_0;...;x^N-x_0]||', '\sigma_{min}^{-1}(L+B)||\epsilon(k)||'}, 'Location', cfg.legend_location);
        else
            lgd1 = legend(ax1, {'||[x^1-x_0;...;x^N-x_0]||'}, 'Location', cfg.legend_location);
        end
        lgd1.FontSize = cfg.legend_font_size;
        lgd1.Interpreter = cfg.interpreter;
        title(ax1, 'Norma de Tracking', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
        xlabel(ax1, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        ylabel(ax1, 'norma 2', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        ax1.FontSize = cfg.axes_font_size;
        ylim(ax1, [0, max(tracking_global)*1.1]);

        ax2 = subplot(1,2,2); hold(ax2, 'on');
        grid(ax2, cfg.grid);
        box(ax2, cfg.box);
        plot(ax2, data.t, eps_global_norm, '-', 'LineWidth', cfg.line_width_closed);
        title(ax2, 'Norma Global de \epsilon(k)', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
        xlabel(ax2, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        ylabel(ax2, '||\epsilon(k)||_2', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        ax2.FontSize = cfg.axes_font_size;
        lgd2 = legend(ax2, {'||[\epsilon_1; ...; \epsilon_N]||'}, 'Location', cfg.legend_location);
        lgd2.FontSize = cfg.legend_font_size;
        lgd2.Interpreter = cfg.interpreter;
    else
        ax = axes(fig);
        hold(ax, 'on');
        grid(ax, cfg.grid);
        box(ax, cfg.box);
        plot(ax, data.t, eps_global_norm, '-', 'LineWidth', cfg.line_width_closed);
        title(ax, 'Norma Global de \epsilon(k)', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
        xlabel(ax, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        ylabel(ax, '||\epsilon(k)||_2', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
        ax.FontSize = cfg.axes_font_size;
        lgd = legend(ax, {'||[\epsilon_1; ...; \epsilon_N]||'}, 'Location', cfg.legend_location);
        lgd.FontSize = cfg.legend_font_size;
        lgd.Interpreter = cfg.interpreter;
    end
end

function fig = plot_control_inputs(data, cfg)
    required_fields = {'t_u', 'u_cl', 'nagent', 'agent_labels'};
    validate_fields(data, required_fields, 'control_inputs');

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);
    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, cfg.grid);
    box(ax, cfg.box);

    for i = 1:data.nagent
        plot(ax, data.t_u, data.u_cl(i,:), '-', 'LineWidth', cfg.line_width_closed);
    end

    title(ax, 'u_i(k) aplicado (malha fechada)', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
    xlabel(ax, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ylabel(ax, 'u_i', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ax.FontSize = cfg.axes_font_size;
    lgd = legend(ax, data.agent_labels, 'Location', cfg.legend_location);
    lgd.FontSize = cfg.legend_font_size;
    lgd.Interpreter = cfg.interpreter;
end

function fig = plot_phase_plane(data, cfg)
    required_fields = {'x_closed', 'nagent', 'agent_labels'};
    validate_fields(data, required_fields, 'phase_plane');
    has_leader = isfield(data, 'leader') && ~isempty(data.leader);

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);
    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, cfg.grid);
    box(ax, cfg.box);

    for i = 1:data.nagent
        plot(ax, squeeze(data.x_closed(1,i,:)), squeeze(data.x_closed(2,i,:)), '-', 'LineWidth', cfg.line_width_closed);
    end
    if has_leader
        plot(ax, data.leader(1,:), data.leader(2,:), 'k:', 'LineWidth', cfg.line_width_reference);
    end

    title(ax, 'Plano de fase: x_1 x x_2 (malha fechada)', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
    xlabel(ax, 'x_1', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ylabel(ax, 'x_2', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ax.FontSize = cfg.axes_font_size;
    legend_labels = get_state_legend_labels(data, has_leader);
    lgd = legend(ax, legend_labels, 'Location', cfg.legend_location);
    lgd.FontSize = cfg.legend_font_size;
    lgd.Interpreter = cfg.interpreter;
end

function fig = plot_phase_3d_time(data, cfg)
    required_fields = {'t', 'x_closed', 'nagent', 'agent_labels'};
    validate_fields(data, required_fields, 'phase_3d_time');
    has_leader = isfield(data, 'leader') && ~isempty(data.leader);

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);
    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, cfg.grid);
    box(ax, cfg.box);

    for i = 1:data.nagent
        plot3(ax, squeeze(data.x_closed(1,i,:)), squeeze(data.x_closed(2,i,:)), data.t(:), '-', 'LineWidth', cfg.line_width_closed);
    end
    if has_leader
        plot3(ax, data.leader(1,:), data.leader(2,:), data.t, 'k:', 'LineWidth', cfg.line_width_reference);
    end

    title(ax, 'Trajetoria 3D: (x_1, x_2, k) - malha fechada', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
    xlabel(ax, 'x_1', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ylabel(ax, 'x_2', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    zlabel(ax, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    view(ax, 45, 25);
    ax.FontSize = cfg.axes_font_size;
    legend_labels = get_state_legend_labels(data, has_leader);
    lgd = legend(ax, legend_labels, 'Location', cfg.legend_location);
    lgd.FontSize = cfg.legend_font_size;
    lgd.Interpreter = cfg.interpreter;
end

function fig = plot_uncertainty_norms(data, cfg)
    required_fields = {'t_u', 'deltaF_norm'};
    validate_fields(data, required_fields, 'uncertainty_norms');

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);
    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, cfg.grid);
    box(ax, cfg.box);

    plot(ax, data.t_u, data.deltaF_norm, 'k-', 'LineWidth', cfg.line_width_reference);

    legends = {'||\DeltaF(k)||'};
    if isfield(data, 'deltaG_norm') && ~isempty(data.deltaG_norm)
        nagent_g = size(data.deltaG_norm, 1);
        for i = 1:nagent_g
            plot(ax, data.t_u, data.deltaG_norm(i,:), '-', 'LineWidth', cfg.line_width_closed);
            legends{end+1} = sprintf('||\\DeltaG_%d(k)||', i); %#ok<AGROW>
        end
    end

    title(ax, 'Norma das incertezas', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
    xlabel(ax, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ylabel(ax, 'norma 2', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ax.FontSize = cfg.axes_font_size;
    lgd = legend(ax, legends, 'Location', cfg.legend_location);
    lgd.FontSize = cfg.legend_font_size;
    lgd.Interpreter = cfg.interpreter;
end

function fig = plot_tracking_global_and_components(data, cfg)
    required_fields = {'t', 'epsilon', 'nstate', 'nagent'};
    validate_fields(data, required_fields, 'tracking_global_and_components');

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);
    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, cfg.grid);
    box(ax, cfg.box);

    comp_labels = cell(1, data.nstate * data.nagent);
    idx_label = 1;
    for i = 1:data.nagent
        for j = 1:data.nstate
            comp_ij = squeeze(data.epsilon(j,i,:));
            plot(ax, data.t, comp_ij, '-', 'LineWidth', cfg.line_width_closed);
            comp_labels{idx_label} = sprintf('\\epsilon_{%d%d}', i, j);
            idx_label = idx_label + 1;
        end
    end

    title(ax, 'Tracking Error Dynamics - Malha Fechada', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
    xlabel(ax, 'Iteration Steps', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ylabel(ax, 'Tracking Error Dynamics', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ax.FontSize = cfg.axes_font_size;
    lgd = legend(ax, comp_labels, 'Location', 'southeast');
    lgd.FontSize = cfg.legend_font_size;
    lgd.Interpreter = cfg.interpreter;
end

function fig = plot_costs(data, cfg)
    required_fields = {'t_u', 'stage_cost_cl', 'cumulative_cost_cl'};
    validate_fields(data, required_fields, 'costs');

    fig = figure('Color', cfg.figure_color, 'Name', cfg.figure_name);

    ax1 = subplot(1,2,1);
    hold(ax1, 'on');
    grid(ax1, cfg.grid);
    box(ax1, cfg.box);
    plot(ax1, data.t_u, data.stage_cost_cl, '-', 'LineWidth', cfg.line_width_closed);
    title(ax1, 'Custo instantaneo', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
    xlabel(ax1, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ylabel(ax1, 'custo', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ax1.FontSize = cfg.axes_font_size;

    ax2 = subplot(1,2,2);
    hold(ax2, 'on');
    grid(ax2, cfg.grid);
    box(ax2, cfg.box);
    plot(ax2, data.t_u, data.cumulative_cost_cl, '-', 'LineWidth', cfg.line_width_closed);
    title(ax2, 'Custo acumulado', 'FontSize', cfg.title_font_size, 'Interpreter', cfg.interpreter);
    xlabel(ax2, 'k', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ylabel(ax2, 'custo', 'FontSize', cfg.label_font_size, 'Interpreter', cfg.interpreter);
    ax2.FontSize = cfg.axes_font_size;
end

function legend_labels = get_state_legend_labels(data, has_leader)
    if isfield(data, 'legend_states') && ~isempty(data.legend_states)
        legend_labels = data.legend_states;
    elseif isfield(data, 'agent_labels') && ~isempty(data.agent_labels)
        legend_labels = data.agent_labels;
    else
        legend_labels = arrayfun(@(i) sprintf('ag%d', i), 1:data.nagent, 'UniformOutput', false);
    end

    if ~has_leader && numel(legend_labels) > data.nagent
        legend_labels = legend_labels(1:data.nagent);
    end
end

function value = get_option(opts, field_name, default_value)
    if isfield(opts, field_name) && ~isempty(opts.(field_name))
        value = opts.(field_name);
    else
        value = default_value;
    end
end

function validate_fields(data, required_fields, plot_name)
    for k = 1:numel(required_fields)
        if ~isfield(data, required_fields{k})
            error('Plot "%s" requer o campo data.%s.', plot_name, required_fields{k});
        end
    end
end
