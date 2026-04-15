clc; clear; %close all;
rng(879797);
%% ======== Config ========
PROBLEM_ID =6;
SAVE_RESULTS = false;
SAVENAME="LQR_Online";
SHOW_ANIMATION = false;   % true/false para ligar/desligar
ANIM_PAUSE_S = 1/144;      % pausa entre frames
RUN_FIXED_GAIN=false;
UPDATE_UNCERTAINTY_EACH_K = false;
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
rest1 = max(abs(1 - (eigL/2)));
Finst= eig(Fk)
Finst=Finst(abs(Finst) > 1);
rest2 = 1/abs(prod(Finst));
fprintf('Restrição de delta: %.4g < delta < %.4g\n', rest1, rest2);
delta=0.35;
%% ======== Parametros RLQR (online) ========
mu = 1e12;
alfa=0.01;
InVal=1;
Q = eye(nstate); 

Q=diag([0.6 0.6]);
Qx=Q;
R = eye(ninput)*InVal;
%calculo de C seguindo o artigo ;

%% ======== inicialização ========
for i = 1:nagent
    P_online{i}=eye(nstate);
end





c=        0.85;
aux_gain=1;


nuncertain = size(Ef,1);
tic()
auxP=[];
for k = 1:300
        for i = 1:nagent
            %Q=delta^2*Fk'*P_online{i}*Gk*((R+Gk'*P_online{i}*Gk)^(-1))*Gk'*P_online{i}*Fk+Qx; 
            [P_online{i}, K_now] = Recursive_RLQR_Norma( ...
                1, mu, alfa, ninput, nstate, nuncertain, ...
                Fk, Gk, Q, R, Hk, Ef, EGk, P_online{i},c);
        end

end
tss=toc();


