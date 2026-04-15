clc; clear; close all;
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

%% ======== Configuracao dos arquivos ========
run_date = '2026-04-08';
file_rlqr = sprintf('resultado_robusto_%s_LQR_Online.mat', run_date);
file_lmi  = sprintf('resultado_robusto_%s_LMI.mat', run_date);
num_iter_plot = 120; % numero maximo de iteracoes a plotar (para evitar sobrecarga em casos muito longos)

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

data_dir = fullfile(script_dir, 'Dados');
path_rlqr = fullfile(data_dir, file_rlqr);
path_lmi  = fullfile(data_dir, file_lmi);

%% ======== Validacao ========
if ~isfile(path_rlqr)
    error('Arquivo RLQR nao encontrado: %s', path_rlqr);
end

if ~isfile(path_lmi)
    error('Arquivo LMI nao encontrado: %s', path_lmi);
end

%% ======== Carga dos resultados ========
data_rlqr = load(path_rlqr, 'results');
data_lmi  = load(path_lmi, 'results');

if ~isfield(data_rlqr, 'results')
    error('O arquivo RLQR nao contem a variavel ''results'': %s', path_rlqr);
end

if ~isfield(data_lmi, 'results')
    error('O arquivo LMI nao contem a variavel ''results'': %s', path_lmi);
end

results_rlqr = data_rlqr.results;
results_lmi  = data_lmi.results;

%% ======== Resumo inicial ========
fprintf('Arquivos carregados com sucesso.\n');
fprintf('RLQR: %s\n', path_rlqr);
fprintf('LMI : %s\n\n', path_lmi);

fprintf('Problema RLQR: %s (ID=%d)\n', results_rlqr.problem_name, results_rlqr.problem_id);
fprintf('Problema LMI : %s (ID=%d)\n', results_lmi.problem_name, results_lmi.problem_id);
fprintf('Horizonte RLQR: T=%d\n', results_rlqr.T);
fprintf('Horizonte LMI : T=%d\n', results_lmi.T);

compare_root_dir = fullfile(script_dir, 'Compare');
if ~exist(compare_root_dir, 'dir')
    mkdir(compare_root_dir);
end

compare_dir = fullfile(compare_root_dir, run_date);
if ~exist(compare_dir, 'dir')
    mkdir(compare_dir);
end

%% ======== Grafico: norma do erro ========
n_iter_error = min([num_iter_plot-90, size(results_rlqr.eps_cl, 3), size(results_lmi.eps_cl, 3)]);
k_error = 0:(n_iter_error-1);

error_norm_rlqr = zeros(1, n_iter_error);
error_norm_lmi = zeros(1, n_iter_error);

for idx = 1:n_iter_error
    eps_rlqr_k = results_rlqr.eps_cl(:, :, idx);
    eps_lmi_k = results_lmi.eps_cl(:, :, idx);
    error_norm_rlqr(idx) = norm(eps_rlqr_k, 'fro');
    error_norm_lmi(idx) = norm(eps_lmi_k, 'fro');
end

fig_error = figure('Color', 'w', 'Position', [100, 100, 1100, 700]);
plot(k_error, error_norm_lmi, 'LineWidth', 6, 'Color', [0.00 0.45 0.74]);
hold on;
plot(k_error, error_norm_rlqr, 'LineWidth',6, 'Color', [0.85 0.33 0.10]);
grid on;
box on;

xlabel('time steps $k$', 'Interpreter', 'latex', 'FontSize', 36, 'FontWeight', 'bold');
ylabel('$\|\varepsilon_k\|$', 'Interpreter', 'latex', 'FontSize', 64, 'FontWeight', 'bold');
legend({'LMI Xu et al 2017. \quad', 'MAS-RLQR \quad'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'northeast');

set(gca, 'FontSize', 36, 'LineWidth', 1.2);
set(fig_error, 'PaperPositionMode', 'auto');

base_name_error = fullfile(compare_dir, 'error_norm_compare');
savefig(fig_error, [base_name_error '.fig']);
exportgraphics(fig_error, [base_name_error '.svg'], 'ContentType', 'vector');
exportgraphics(fig_error, [base_name_error '.pdf'], 'ContentType', 'vector');

fprintf('\nFigura salva em:\n');
fprintf('%s.fig\n', base_name_error);
fprintf('%s.svg\n', base_name_error);
fprintf('%s.pdf\n', base_name_error);

%% ======== Grafico: custo acumulado ========
n_iter_cost = min([num_iter_plot-70, numel(results_rlqr.cumulative_cost_cl), numel(results_lmi.cumulative_cost_cl)]);
k_cost = 0:(n_iter_cost-1);

cost_rlqr = results_rlqr.cumulative_cost_cl(1:n_iter_cost);
cost_lmi = results_lmi.cumulative_cost_cl(1:n_iter_cost);

fig_cost = figure('Color', 'w', 'Position', [150, 150, 1100, 700]);
plot(k_cost, cost_lmi, 'LineWidth', 4, 'Color', [0.00 0.45 0.74]);
hold on;
plot(k_cost, cost_rlqr, 'LineWidth',4, 'Color', [0.85 0.33 0.10]);
grid on;
box on;

xlabel('time steps $k$', 'Interpreter', 'latex', 'FontSize', 36, 'FontWeight', 'bold');
ylabel('$J_k$', 'Interpreter', 'latex', 'FontSize', 64, 'FontWeight', 'bold');
legend({'LMI Xu et al 2017. \quad', 'MAS-RLQR \quad'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'southeast');

set(gca, 'FontSize', 36, 'LineWidth', 1.2);
set(fig_cost, 'PaperPositionMode', 'auto');

base_name_cost = fullfile(compare_dir, 'cumulative_cost_compare');
savefig(fig_cost, [base_name_cost '.fig']);
exportgraphics(fig_cost, [base_name_cost '.svg'], 'ContentType', 'vector');
exportgraphics(fig_cost, [base_name_cost '.pdf'], 'ContentType', 'vector');

fprintf('\nFigura salva em:\n');
fprintf('%s.fig\n', base_name_cost);
fprintf('%s.svg\n', base_name_cost);
fprintf('%s.pdf\n', base_name_cost);

%% ======== Graficos: trajetorias RLQR por estado ========
n_iter_state = min([num_iter_plot, size(results_rlqr.x_cl, 3)]);
k_state = 0:(n_iter_state-1);
num_agents = size(results_rlqr.x_cl, 2);
agent_colors = lines(num_agents);

fig_x1 = figure('Color', 'w', 'Position', [200, 200, 1100, 700]);
hold on;
for agent_idx = 1:num_agents
    x1_agent = squeeze(results_rlqr.x_cl(1, agent_idx, 1:n_iter_state));
    plot(k_state, x1_agent, 'LineWidth',4, 'Color', agent_colors(agent_idx, :));
end
grid on;
box on;

xlabel('time steps $k$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('$x_{1,k}^i$', 'Interpreter', 'latex', 'FontSize', 48, 'FontWeight', 'bold');
legend(arrayfun(@(i) sprintf('Agent %d', i), 1:num_agents, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 24, 'Location', 'best');

set(gca, 'FontSize', 36, 'LineWidth', 1.2);
set(fig_x1, 'PaperPositionMode', 'auto');

base_name_x1 = fullfile(compare_dir, 'rlqr_state1_agents');
savefig(fig_x1, [base_name_x1 '.fig']);
exportgraphics(fig_x1, [base_name_x1 '.svg'], 'ContentType', 'vector');
exportgraphics(fig_x1, [base_name_x1 '.pdf'], 'ContentType', 'vector');

fprintf('\nFigura salva em:\n');
fprintf('%s.fig\n', base_name_x1);
fprintf('%s.svg\n', base_name_x1);
fprintf('%s.pdf\n', base_name_x1);

fig_x2 = figure('Color', 'w', 'Position', [250, 250, 1100, 700]);
hold on;
for agent_idx = 1:num_agents
    x2_agent = squeeze(results_rlqr.x_cl(2, agent_idx, 1:n_iter_state));
    plot(k_state, x2_agent, 'LineWidth', 4, 'Color', agent_colors(agent_idx, :));
end
grid on;
box on;

xlabel('time steps $k$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('$x_{2,k}^i$', 'Interpreter', 'latex', 'FontSize', 48, 'FontWeight', 'bold');
legend(arrayfun(@(i) sprintf('Agent %d', i), 1:num_agents, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 24, 'Location', 'best');

set(gca, 'FontSize', 36, 'LineWidth', 1.2);
set(fig_x2, 'PaperPositionMode', 'auto');

base_name_x2 = fullfile(compare_dir, 'rlqr_state2_agents');
savefig(fig_x2, [base_name_x2 '.fig']);
exportgraphics(fig_x2, [base_name_x2 '.svg'], 'ContentType', 'vector');
exportgraphics(fig_x2, [base_name_x2 '.pdf'], 'ContentType', 'vector');

fprintf('\nFigura salva em:\n');
fprintf('%s.fig\n', base_name_x2);
fprintf('%s.svg\n', base_name_x2);
fprintf('%s.pdf\n', base_name_x2);
