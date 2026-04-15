clc; clear; close all;
rng(879797);


SAVE_RESULTS = false;

%% ======== Config ========
PROBLEM_ID = 1;
SAVE_RESULTS = false;
SAVENAME="Exemplo1_robustoComparar";
SHOW_ANIMATION = false;   % true/false para ligar/desligar
ANIM_PAUSE_S = 1/144;      % pausa entre frames
RUN_FIXED_GAIN=false;
UPDATE_UNCERTAINTY_EACH_K = true;
LQRONLINE=true;


%% ======== Dados via banco unificado ========

problems = Consensus_Problems();

if PROBLEM_ID > numel(problems) || isempty(problems{PROBLEM_ID})
    error('PROBLEM_ID=%d nao esta definido em Consensus_Problems.m', PROBLEM_ID);
end

p = problems{PROBLEM_ID};
A = p.F;
B=p.Gk;
Lcal = p.Lcal;
E1=p.Ef;
E2=p.EGk;
D=p.H;


nagent = size(Lcal,1);
nstate = size(A,1);
ninput = size(B,2);
lam = eig(Lcal); lam = lam(abs(lam) > 1e-12);
Qx=diag([1 1]);
Qu= eye(ninput);
fprintf('lambda_max(Î“) = %.4f\n', max(real(lam)));
tol = 0;
%% parametros para LMI
EigVal=eig(Lcal);
EigVal=EigVal(abs(EigVal)>1e-12);

    Ai = blkdiag(A, A);
    Bi = blkdiag(B, B);
   
    Di = blkdiag(D, D);
    E1i = blkdiag(E1, E1);
    E2i = blkdiag(E2, E2);
    Qxi = blkdiag(Qx, Qx);
    Qui = blkdiag(Qu, Qu);
    
 
    Quiinv=inv(Qui);
    epsi = sdpvar(1,1,'full');
    Xi=sdpvar(nstate*2,nstate*2,'sym');
    Wi=sdpvar(ninput*2,nstate*2,'full');
    LMIs=[epsi>=1e-6, Xi>=1e-6*eye(2*nstate)];
for i = 1:size(EigVal,1)
    ri{i} = [real(EigVal(i))*eye(nstate) -imag(EigVal(i)*eye(nstate))
    imag(EigVal(i))*eye(nstate) real(EigVal(i))*eye(nstate) ];
    gammai{i} = [real(EigVal(i))*eye(ninput) -imag(EigVal(i)*eye(ninput))
    imag(EigVal(i))*eye(ninput) real(EigVal(i))*eye(ninput)];
    Qx_i = ri{i}'*Qxi*ri{i};
Qx_i = (Qx_i + Qx_i')/2;
Qxiinv = inv(Qx_i);
    nx = size(Xi,1);      % 2d
    nu = size(Qui,1);     % 2m
    nw = size(Di,2);      % dimensão do F
    Ieps = eye(nw);
    
    M11 = epsi*Di*Di' - Xi;
    M12 = Ai*Xi - ri{i}*Bi*Wi;
    M13 = zeros(nx,nw);
    M14 = zeros(nx,nx);
    M15 = zeros(nx,nu);

    M22 = -Xi;
    M23 = (E1i*Xi - ri{i}*E2i*Wi)';
    M24 = Xi;
    M25 = Wi'*gammai{i}';

    M33 = -epsi*Ieps;
    M34 = zeros(nw,nx);
    M35 = zeros(nw,nu);

    M44 = -Qxiinv;
    M45 = zeros(nx,nu);

    M55 = -Quiinv;

    BigLMI = [M11   M12   M13   M14   M15;
          M12'  M22   M23   M24   M25;
          M13'  M23'  M33   M34   M35;
          M14'  M24'  M34'  M44   M45;
          M15'  M25'  M35'  M45'  M55];
BigLMI = (BigLMI + BigLMI')/2;
BigLMIs{i} =BigLMI;
LMIs = [LMIs, BigLMI <= -1e-8*eye(size(BigLMI,1))];
    end
%% solve LMIs
tic;
ops = sdpsettings('solver','sedumi','verbose',1);
sol = optimize(LMIs, [], ops);
t=toc;
disp(sol.problem)
disp(sol.info)

%% obter K 
Xi_val = value(Xi);
    Wi_val = value(Wi);
     Ktil = Wi_val / Xi_val;
      m = ninput;
    d = nstate;

    K11 = Ktil(1:m,     1:d);
    K22 = Ktil(m+1:2*m, d+1:2*d);
        Kf = 0.5*(K11 + K22);
