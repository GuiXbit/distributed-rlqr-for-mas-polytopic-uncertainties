clc; clear; %close all;
rng(879797);
%% ======== Config ========
PROBLEM_ID =1;
SAVE_RESULTS = false;
SAVENAME="RLQR";
SHOW_ANIMATION = false;   % true/false para ligar/desligar
ANIM_PAUSE_S = 1/144;      % pausa entre frames
RUN_FIXED_GAIN=false;
UPDATE_UNCERTAINTY_EACH_K = true;
PLOTS=true;

   
%Applied Control
LQRONLINE=true;
LQROFFLINE=false;
LMI=false;
consenso2=false;

%% ======== Dados via banco unificado ========
problems = Consensus_Problems();

if PROBLEM_ID > numel(problems) || isempty(problems{PROBLEM_ID})
    error('PROBLEM_ID=%d nao esta definido em Consensus_Problems.m', PROBLEM_ID);
end

p = problems{PROBLEM_ID};
Fk = p.F;
Gk=p.Gk;
Lcal = p.Lcal;
Ef=p.Ef;
EGk=p.EGk;
Hk=p.H;


nagent = size(Lcal,1);
nstate = size(Fk,1);
ninput = size(Gk,2);

auxlaplace=[];
for i=1:nagent
    auxlaplace=[auxlaplace;Lcal(i,:)/Lcal(i,i)];
end
lam = eig(auxlaplace); lam = lam(abs(lam) > 1e-12);
fprintf('lambda_max(Î“) = %.4f\n', max(real(lam)));
eigL = lam;
rest1 = max(abs(1 - (eigL)));
Finst= eig(Fk)
Finst=Finst(abs(Finst) > 1);
rest2 = 1/abs(prod(Finst));
fprintf('Restrição de delta: %.4g < delta < %.4g\n', rest1, rest2);
delta=0.9;
%% ======== Parametros RLQR (online) ========
mu = 1e12;
alfa=0.01;
InVal=1;
Q = eye(nstate); 

Q=diag([1 1]);
Qx=Q;
R = eye(ninput)*InVal;
%calculo de C seguindo o artigo ;

%% ======== inicialização ========
T = p.T;
% estados Riccati para atualizacao online
P_online = cell(nagent,1);
K_last = cell(nagent,1);
for i = 1:nagent
    P_online{i} = eye(nstate);
    K_last{i} = zeros(size(Gk,2), nstate);
end

s = size(Hk, 2);
l1 = size(Ef, 1);
l2 = size(EGk, 1);
normalize_uncertainty = @(M) M / max(1, norm(M, 2));

eigLNS=eig(Lcal); 
eigLNS=eigLNS(abs(eigLNS)>1e-12);
r1=max(eigLNS)/min(eigLNS);

% condições iniciais dos agentes
if ~isfield(p,'x_init') || isempty(p.x_init)
    error('O problema %d nao possui campo x_init em Consensus_Problems.m', PROBLEM_ID);
else
    x_init = p.x_init;
    if size(x_init,1) ~= nstate || size(x_init,2) ~= nagent
        error('x_init deve ter dimensao %dx%d para o problema %d', nstate, nagent, PROBLEM_ID);
    end
end

x_open  = zeros(nstate,nagent,T+1);
x_cl    = zeros(nstate,nagent,T+1);

eps_open = zeros(nstate,nagent,T+1);
eps_cl   = zeros(nstate,nagent,T+1);
deltaF_norm = zeros(1,T);
deltaG_norm = zeros(1,T);

u_open = zeros(nagent,T);
u_cl   = zeros(nagent,T);
stage_cost_cl = zeros(1,T);
cumulative_cost_cl = zeros(1,T);
stage_cost_cl_agents = zeros(nagent,T);

x_open(:,:,1) = x_init;
x_cl(:,:,1)   = x_init;

% ---- calcula Y_i = sum_j a_ij (x_j - x_i) ----
calc_Yk = @(X) ...
    reshape(-kron(Lcal, eye(nstate)) * reshape( X, [], 1 ), nstate, nagent);
lam = eig(Lcal);
lam = lam(abs(lam) > 1e-12);

maxeigLNS=max(abs(lam));
mineigLNS=min(abs(eigLNS));
c=       1;
aux_gain=1;
lambda=(1 + alfa) * norm( mu * Hk' * Hk );
aux=Gk*(R^(-1))*EGk'*(lambda^(-1)*eye(size(EGk,1))+EGk*(R^(-1))*EGk')^(-1)*Ef;
Calf=Fk-aux;
Finst= eig(Calf)
Finst=Finst(abs(Finst) > 1);
rest2 = 1/abs(prod(Finst));
fprintf('Restrição de delta: %.4g < delta < %.4g\n', rest1, rest2);
%% ======== Loop de simulaÃ§Ã£o ========
if ~UPDATE_UNCERTAINTY_EACH_K
  DELTA=diag([0.35 0.25 0.15]);
    deltaF = Hk * DELTA * Ef;
    deltaG = Hk * DELTA * EGk;
    DeltaG_base=DELTA;
    DELTAF= DELTA;
end
tic()
for k = 1:T
    if (UPDATE_UNCERTAINTY_EACH_K && (mod(k,1)==0)|| k==1 && UPDATE_UNCERTAINTY_EACH_K)
        DELTAF = CalcDelta(s, l1,"aaa",k)*0.8;  % Random matrix in [-0.5, 0.5]
        DELTAF = normalize_uncertainty(DELTAF);
        deltaF = Hk * DELTAF * Ef;         
        DeltaG_base = CalcDelta(s, l2,"sinrand",k)*0;  % Random matrix in [-1, 1]
        DeltaG_base = normalize_uncertainty(DeltaG_base);
        DeltaG_base=DELTAF;
        deltaG = Hk * DeltaG_base * EGk;
        
    end
    deltaF_norm(k) = norm(DELTAF, 2);
    deltaG_norm(k) = norm(DeltaG_base, 2);
    % --- malha aberta ---
    E = calc_Yk(x_open(:,:,k));
    eps_open(:,:,k) = E;
    for i = 1:nagent
        x_open(:,i,k+1) = (Fk+deltaF)*x_open(:,i,k);
    end
    % --- malha fechada (controle via epsilon) ---
    E = calc_Yk(x_cl(:,:,k));
    eps_cl(:,:,k) = E;
    nuncertain = size(Ef,1);

    
        for i = 1:nagent
            if LQRONLINE
            Q=delta^2*Fk'*P_online{i}*Gk*((R+Gk'*P_online{i}*Gk)^(-1))*Gk'*P_online{i}*Fk+Qx; 
            [P_online{i}, K_now] = Recursive_RLQR_Pol( ...
                1, mu, alfa, ninput, nstate, nuncertain, ...
                Fk, Gk, Q, R, Hk, Ef, EGk, P_online{i},c);
                aux_gain=Lcal(i,i);
            elseif LMI
                 K_now=[   0.2132    0.5146];
            elseif LQROFFLINE
                K_now=[  0.0585   -0.2499 0.6732];
                aux_gain=Lcal(i,i);
            elseif consenso2
                [P_online{i}, K_now] = Recursive_Lqr_Consensus(Fk, Gk, Q, R, P_online{i},delta);
                aux_gain=Lcal(i,i);
            end
            K_last{i} = K_now;             % ou K_now(1:ninput,:), se precisar
            u_i = K_now * E(:,i);       % ou -Kself * E(:,i), dependendo da sua funcao
            u_cl(i,k) = norm(u_i,2);
            x_cl(:,i,k+1) = (Fk + deltaF) * x_cl(:,i,k) + (Gk + deltaG) * u_i/aux_gain;
            stage_cost_cl_agents(i,k) = E(:,i)' * Qx * E(:,i) + u_i' * R * u_i/aux_gain;
        end
        stage_cost_cl(k) = sum(stage_cost_cl_agents(:,k));
        cumulative_cost_cl(k) = sum(stage_cost_cl(1:k));
end
tss=toc();
for i = 1:nagent
    eigVals = eig((Fk + deltaF) + (Gk + deltaG) * K_last{i});
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
PLOT_ZOOM.global_tracking.enable = false;
PLOT_ZOOM.global_tracking.xlim = [36 38];
PLOT_ZOOM.global_tracking.ylim = [0.5 0.75];
PLOT_ZOOM.global_tracking.position = [0.20 0.30 0.18 0.18];
PLOT_ZOOM.global_eta.enable = false;
PLOT_ZOOM.global_eta.xlim = [30 50];
PLOT_ZOOM.global_eta.ylim = [];
PLOT_ZOOM.global_eta.position = [0.58 0.52 0.30 0.30];

% 1) TrajetÃ³ria dos estados (ex.: componente 1)
% Plot 1: trajetorias dos estados em malha aberta e fechada.
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
plotData.deltaF_norm = deltaF_norm;
plotData.deltaG_norm = deltaG_norm;
plotData.stage_cost_cl = stage_cost_cl;
plotData.cumulative_cost_cl = cumulative_cost_cl;


if PLOTS
    plotOptsStates = struct('figure_name', 'Trajetorias de estados (aberta x fechada)');
    PlotsLQR_Reform('state_trajectories', plotData, plotOptsStates);

    plotOptsTracking = struct('figure_name', 'Tracking error por componentes');
    PlotsLQR_Reform('epsilon_components', plotData, plotOptsTracking);

    plotOptsEpsNorm = struct('figure_name', 'Norma global de epsilon');
    PlotsLQR_Reform('epsilon_global_norm', plotData, plotOptsEpsNorm);

    plotOptsControl = struct('figure_name', 'Controles aplicados');
    PlotsLQR_Reform('control_inputs', plotData, plotOptsControl);

    plotOptsUnc = struct('figure_name', 'Norma das incertezas');
    PlotsLQR_Reform('uncertainty_norms', plotData, plotOptsUnc);

    plotOptsPhase = struct('figure_name', 'Plano de fase');
    PlotsLQR_Reform('phase_plane', plotData, plotOptsPhase);

    plotOptsPhase3D = struct('figure_name', 'Trajetoria 3D');
    PlotsLQR_Reform('phase_3d_time', plotData, plotOptsPhase3D);

    plotOptsCost = struct('figure_name', 'Custos (instantaneo e acumulado)');
    PlotsLQR_Reform('costs', plotData, plotOptsCost);
end

%% 6) Animacao 

if SHOW_ANIMATION
    figAnim = figure('Name','Animacao consenso (x1x2, 3D, erro)');
    tl = tiledlayout(figAnim,1,3,'Padding','compact','TileSpacing','compact');

    % --- painel 1: plano de fase x1 x2 ---
    ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
    title(ax1,'Plano de fase x_1 x x_2');
    xlabel(ax1,'x_1'); ylabel(ax1,'x_2');
    hPhase = gobjects(nagent,1);
    for i = 1:nagent
        hPhase(i) = plot(ax1, nan, nan, '-', 'LineWidth', 1.4);
    end

    % --- painel 2: 3D (x1, x2, k) ---
    ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
    title(ax2,'3D: (x_1, x_2, k)');
    xlabel(ax2,'x_1'); ylabel(ax2,'x_2'); zlabel(ax2,'k');
    view(ax2,45,25);
    h3D = gobjects(nagent,1);
    for i = 1:nagent
        h3D(i) = plot3(ax2, nan, nan, nan, '-', 'LineWidth', 1.3);
    end

    % --- painel 3: norma do erro ||epsilon_i(k)|| ---
    ax3 = nexttile(tl,3); hold(ax3,'on'); grid(ax3,'on');
    title(ax3,'Norma do erro ||\epsilon_i(k)||');
    xlabel(ax3,'k'); ylabel(ax3,'norma 2');
    hErr = gobjects(nagent,1);
    for i = 1:nagent
        hErr(i) = plot(ax3, nan, nan, '-', 'LineWidth', 1.3);
    end
    legend(ax3, agentLabels, 'Location', 'best');

    % limites fixos para evitar autoscale durante a animacao
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
    output_file = fullfile(output_dir, sprintf('resultado_robusto_%s_%s.mat', run_date, SAVENAME));

    results = struct();
    results.mode = 'robusto';
    results.problem_id = PROBLEM_ID;
    results.problem_name = p.name;
    results.run_date = run_date;
    results.rng_seed = 879797;
    results.Fk = Fk;
    results.Lcal = Lcal;
    results.Gk = Gk;
    results.Hk = Hk;
    results.EGk = EGk;
    results.Ef = Ef;
    results.T = T;
    results.x_init = x_init;
    results.mu = mu;
    results.alfa = alfa;
    results.Q = Q;
    results.x_open = x_open;
    results.x_cl = x_cl;
    results.eps_open = eps_open;
    results.eps_cl = eps_cl;
    results.u_open = u_open;
    results.u_cl = u_cl;
    results.plot_data = plotData;
    results.stage_cost_cl = stage_cost_cl;
    results.cumulative_cost_cl = cumulative_cost_cl;
    results.stage_cost_cl_agents = stage_cost_cl_agents;
    results.P_online = P_online;
    results.K_last = K_last;
    results.deltaF = deltaF;
    results.deltaG = deltaG;
    results.DELTAF = DELTAF;
    results.DeltaG_base = DeltaG_base;
    results.deltaF_norm = deltaF_norm;
    results.deltaG_norm = deltaG_norm;
    results.UPDATE_UNCERTAINTY_EACH_K = UPDATE_UNCERTAINTY_EACH_K;

    save(output_file, 'results');
    fprintf('Resultados salvos em: %s\n', output_file);
end



% % Plot 5: Custo
% figure;
% subplot(1,2,1); hold on; grid on;
% plot(0:T-1, stage_cost_cl, '-', 'LineWidth', 1.8);
% if RUN_FIXED_GAIN
%     plot(0:T-1, stage_cost_cl_fixed, '--', 'LineWidth', 1.8);
% end
% title('Custo instantaneo');
% xlabel('k');
% ylabel('custo');
% if RUN_FIXED_GAIN
%     legend({'online', 'K_{lmi} fixo'}, 'Location', 'best');
% else
%     legend({'online'}, 'Location', 'best');
% end

% subplot(1,2,2); hold on; grid on;
% plot(0:T-1, cumulative_cost_cl, '-', 'LineWidth', 1.8);
% if RUN_FIXED_GAIN
%     plot(0:T-1, cumulative_cost_cl_fixed, '--', 'LineWidth', 1.8);
% end
% title('Custo acumulado');
% xlabel('k');
% ylabel('custo');
% if RUN_FIXED_GAIN
%     legend({'online', 'K_{lmi} fixo'}, 'Location', 'best');
% else
%     legend({'online'}, 'Location', 'best');
% end


%%
P_online{i}
M=cell(nagent,1);
for i=1:nagent-1
   M{i}=(Fk)'*P_online{i}*Gk*(R+Gk'*P_online{i}*Gk)^(-1)*Gk'*P_online{i}*Fk;
   C{i}=(1-1/max(eig((sqrt(Q)^(-1))*M{i}*(sqrt(Q)^(-1)))));
   X{i}=roots([(eigLNS(i)^2)/4 -real(eigLNS(i)) C{i}]);
   Md{i}=(Fk+deltaF)'*P_online{i}*(Gk+deltaG)*(R+(Gk+deltaG)'*P_online{i}*(Gk+deltaG))^(-1)*(Gk+deltaG)'*P_online{i}*(Fk+deltaF);
   Cd{i}=(1-1/max(eig((sqrt(Q)^(-1))*Md{i}*(sqrt(Q)^(-1)))));
   Xd{i}=roots([(eigLNS(i)^2)/4 -real(eigLNS(i)) Cd{i}]);
   
end