clc; clear; close all;

%% ======== Dados via banco unificado ========
PROBLEM_ID = 2;
SAVE_RESULTS = false;
problems = Consensus_Problems();

if PROBLEM_ID > numel(problems) || isempty(problems{PROBLEM_ID})
    error('PROBLEM_ID=%d nao esta definido em Consensus_Problems.m', PROBLEM_ID);
end

p = problems{PROBLEM_ID};
Fk = p.F;
Lcal = p.Lcal;
Gk = p.Gk;

nagent = size(Lcal,1);
nstate = size(Fk,1);
ninput = size(Gk,2);

%% ======== Parametros RLQR (online) ========
mu = 1e12;
Q = eye(nstate);
R = eye(ninput)*10^4;

%% ======== Simulacao ========
T = p.T;
P_online = cell(nagent,1);
K_last = cell(nagent,1);
for i = 1:nagent
    P_online{i} = eye(nstate);
    K_last{i} = zeros(ninput, nstate);
end

if ~isfield(p,'x_init') || isempty(p.x_init)
    error('O problema %d nao possui campo x_init em Consensus_Problems.m', PROBLEM_ID);
else
    x_init = p.x_init;
    if size(x_init,1) ~= nstate || size(x_init,2) ~= nagent
        error('x_init deve ter dimensao %dx%d para o problema %d', nstate, nagent, PROBLEM_ID);
    end
end
eigLNS=eig(Lcal); 
eigLNS=eigLNS(abs(eigLNS)>1e-12);
lambda_max = max(abs(eigLNS));
lambda_min = min(abs(eigLNS));
r1=lambda_max/lambda_min;


x_open  = zeros(nstate,nagent,T+1);
x_cl    = zeros(nstate,nagent,T+1);

eps_open = zeros(nstate,nagent,T+1);
eps_cl   = zeros(nstate,nagent,T+1);

u_open = zeros(nagent,T);
u_cl   = zeros(nagent,T);

x_open(:,:,1) = x_init;
x_cl(:,:,1)   = x_init;

calc_Yk = @(X) ...
    reshape(-kron(Lcal, eye(nstate)) * reshape(X, [], 1), nstate, nagent);

%% ======== Loop de simulacao ========
for k = 1:T
    E = calc_Yk(x_open(:,:,k));
    eps_open(:,:,k) = E;

    for i = 1:nagent
        x_open(:,i,k+1) = Fk * x_open(:,i,k);
    end

    E = calc_Yk(x_cl(:,:,k));
    eps_cl(:,:,k) = E;

    for i = 1:nagent
        [P_online{i}, K_now, L_now] = Recursive_LQR_c(1, mu, ninput, nstate, Fk, Gk, Q, R, P_online{i},2.7);
        K_last{i} = K_now;
        u_i = K_now * E(:,i);
        u_cl(i,k) = norm(u_i, 2);
        x_cl(:,i,k+1) = Fk * x_cl(:,i,k) + Gk * u_i;
        if k == 1
            eps_clMF(:,:,k) = E;
        end
        rtheta=sqrt(10^4/(1e4+max(eig(Gk'*P_online{i}*Gk))));
        r2=1/(1-rtheta^2);
        eps_clMF(:,i,k+1) = L_now * eps_clMF(:,i,k);
    end
    if  k<10 || mod(k,T/10)==0
        % fprintf(' %s < %s\n', r1, r2);
        % fprintf(' %s < C <%s\n', (1-rtheta)/lambda_min, (1+rtheta)/lambda_max);
    end
end

for i = 1:nagent
    eigVals = eig(Fk + Gk * K_last{i});
    disp(['i = ', num2str(i), '  eig(F + G*K) online final:']);
    disp(eigVals);
    fprintf('  raio espectral (online final) = %.4f\n', max(abs(eigVals)));
end

eps_open(:,:,T+1) = calc_Yk(x_open(:,:,T+1));
eps_cl(:,:,T+1)   = calc_Yk(x_cl(:,:,T+1));

%% ======== Plots ========

t = 0:T;
agentLabels = arrayfun(@(i) sprintf('ag%d', i), 1:nagent, 'UniformOutput', false);
legendStates = agentLabels;

plotData = struct();
plotData.t = t;
plotData.t_u = 0:T-1;
plotData.x_open = x_open;
plotData.x_closed = x_cl;
plotData.epsilon = eps_cl;
plotData.u_cl = u_cl;
plotData.nstate = nstate;
plotData.nagent = nagent;
plotData.agent_labels = agentLabels;
plotData.legend_states = legendStates;

plotOptsStates = struct('figure_name', 'Trajetorias de estados (aberta x fechada)');
PlotsLQR_Reform('state_trajectories', plotData, plotOptsStates);

plotOptsTracking = struct('figure_name', 'Tracking error por componentes');
PlotsLQR_Reform('epsilon_components', plotData, plotOptsTracking);

plotOptsEpsNorm = struct('figure_name', 'Norma global de epsilon');
PlotsLQR_Reform('epsilon_global_norm', plotData, plotOptsEpsNorm);

plotOptsControl = struct('figure_name', 'Controles aplicados');
PlotsLQR_Reform('control_inputs', plotData, plotOptsControl);

plotOptsPhase = struct('figure_name', 'Plano de fase');
PlotsLQR_Reform('phase_plane', plotData, plotOptsPhase);

plotOptsPhase3D = struct('figure_name', 'Trajetoria 3D');
PlotsLQR_Reform('phase_3d_time', plotData, plotOptsPhase3D);

%% 6) Animacao (opcional): fase, 3D e norma do erro evoluindo no tempo
SHOW_ANIMATION = false;
ANIM_PAUSE_S = 1/144;

if SHOW_ANIMATION
    figAnim = figure('Name','Animacao consenso (x1x2, 3D, erro)');
    tl = tiledlayout(figAnim,1,3,'Padding','compact','TileSpacing','compact');

    ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
    title(ax1,'Plano de fase x_1 x x_2');
    xlabel(ax1,'x_1'); ylabel(ax1,'x_2');
    hPhase = gobjects(nagent,1);
    for i = 1:nagent
        hPhase(i) = plot(ax1, nan, nan, '-', 'LineWidth', 1.4);
    end

    ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
    title(ax2,'3D: (x_1, x_2, k)');
    xlabel(ax2,'x_1'); ylabel(ax2,'x_2'); zlabel(ax2,'k');
    view(ax2,45,25);
    h3D = gobjects(nagent,1);
    for i = 1:nagent
        h3D(i) = plot3(ax2, nan, nan, nan, '-', 'LineWidth', 1.3);
    end

    ax3 = nexttile(tl,3); hold(ax3,'on'); grid(ax3,'on');
    title(ax3,'Norma do erro ||\epsilon_i(k)||');
    xlabel(ax3,'k'); ylabel(ax3,'norma 2');
    hErr = gobjects(nagent,1);
    for i = 1:nagent
        hErr(i) = plot(ax3, nan, nan, '-', 'LineWidth', 1.3);
    end
    legend(ax3, agentLabels, 'Location', 'best');

    x1_all = reshape(x_cl(1,:,:),1,[]);
    x2_all = reshape(x_cl(2,:,:),1,[]);
    e_all = zeros(nagent, numel(t));
    for i = 1:nagent
        e_all(i,:) = vecnorm(squeeze(eps_cl(:,i,:)), 2, 1);
    end
    e_max = max(e_all, [], 'all');
    if ~isfinite(e_max) || e_max <= 0
        e_max = 1;
    end
    dx1 = max(1e-6, 0.05*(max(x1_all)-min(x1_all)));
    dx2 = max(1e-6, 0.05*(max(x2_all)-min(x2_all)));
    xlim(ax1, [min(x1_all)-dx1, max(x1_all)+dx1]);
    ylim(ax1, [min(x2_all)-dx2, max(x2_all)+dx2]);
    xlim(ax2, [min(x1_all)-dx1, max(x1_all)+dx1]);
    ylim(ax2, [min(x2_all)-dx2, max(x2_all)+dx2]);
    zlim(ax2, [t(1), t(end)]);
    xlim(ax3, [t(1), t(end)]);
    ylim(ax3, [0, 1.05*e_max]);

    for kk = 1:numel(t)
        for i = 1:nagent
            set(hPhase(i), 'XData', squeeze(x_cl(1,i,1:kk)), 'YData', squeeze(x_cl(2,i,1:kk)));
            set(h3D(i), 'XData', squeeze(x_cl(1,i,1:kk)), ...
                        'YData', squeeze(x_cl(2,i,1:kk)), ...
                        'ZData', t(1:kk));
            e_hist = vecnorm(squeeze(eps_cl(:,i,1:kk)), 2, 1);
            set(hErr(i), 'XData', t(1:kk), 'YData', e_hist);
        end

        title(ax1, sprintf('Plano de fase x_1 x x_2 (k = %d)', t(kk)));
        title(ax2, sprintf('3D: (x_1, x_2, k) (k = %d)', t(kk)));
        title(ax3, sprintf('Norma do erro ||\\epsilon_i(k)|| (k = %d)', t(kk)));

        drawnow;
        pause(ANIM_PAUSE_S);
    end
end

%% ======== Salvamento dos resultados ========
if SAVE_RESULTS
    script_dir = fileparts(mfilename('fullpath'));
    if isempty(script_dir)
        script_dir = pwd;
    end

    output_dir = fullfile(script_dir, 'Dados');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    run_date = datestr(now, 'yyyy-mm-dd');
    output_file = fullfile(output_dir, sprintf('resultado_%s_problem_%02d.mat', run_date, PROBLEM_ID));

    results = struct();
    results.mode = 'nominal';
    results.problem_id = PROBLEM_ID;
    results.problem_name = p.name;
    results.run_date = run_date;
    results.Fk = Fk;
    results.Lcal = Lcal;
    results.Gk = Gk;
    results.T = T;
    results.x_init = x_init;
    results.mu = mu;
    results.Q = Q;
    results.R = R;
    results.x_open = x_open;
    results.x_cl = x_cl;
    results.eps_open = eps_open;
    results.eps_cl = eps_cl;
    results.u_open = u_open;
    results.u_cl = u_cl;
    results.plot_data = plotData;
    results.P_online = P_online;
    results.K_last = K_last;

    save(output_file, 'results');
    fprintf('Resultados salvos em: %s\n', output_file);
end
