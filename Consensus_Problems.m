function problems = Consensus_Problems()
% Banco unificado: cada entrada inclui dinamica + topologia + iniciais.
% Campos usados no script:
%   name, A, Lcal, Bcal, G0 ou Gk, x_init, x0_leader, T

problems = cell(1,8);

% ===== ID 1:  =====
p = struct();
p.name = 'Exemplo 4 agentes artigo robusto';
p.F = [1 1;
       0 1];
p.Ef = [0.2 0.1 ;
        0.3 0.4];
p.Lcal = [ 2  -1 -1 0;
          -1 1 0 0;
           0 -1  2 -1;
           -1 0 -1 2];
p.Gk = [0.5 1]';
p.EGk = [0.2;
         0.4];
p.H=[0.2 0; 0 0.3]';
p.x_init = [ -2  4 -6 -3;
            3  -5  2 -2];
p.T = 140;
problems{1} = p;

% ===== ID 2:  =====
p = struct();
p.name = 'Exemplo 4 agentes artigo C Gain';
p.F = [1 1 0 ;
       0 1 0 ;
       0 0 1];
p.Gk = [0 0
1 0 
0 3];

p.Lcal = [ 0 0 0 0;
          -1 1 0 0;
           0 -1  2 -1;
           0 0 -1 1];

       
p.Ef = [0.1 0 0 ;
        0 0.1 0 ;
        0 0 0.1];
p.EGk = [0.1 0  ;
        0 0.1  ;
        0  0.1];
p.H=[1 0 1 
0 1 0 
1 0 1];



p.x_init = rand(3,4);
p.T = 200;
problems{2} = p;


% ===== ID 3:  =====
p = struct();
p.name = 'Exemplo 4 agentes artigo C Gain';
p.F = [0.7 -0.7 0;
0.7 0.7 0;
0 0 1];
p.Gk = [0.2 -0.4 1]';

p.Lcal = [ 0.7 0 0 0 0 0 0 -0.7
           -1.3 1.3 0 0 0 0 0 0
           0 -1.2 2.3 -1.1 0 0 0 0
           0 0 0 0 0 0 0 0
           0 0 0  -1.2 1.2 0 0 0
           0 0 0 -0.8 -0.9 1.7 0 0
           -0.7 0 0 0 0 -1.4 2.1 0
          0 0 0 0 0 0 -1 1];
       
p.Ef = [0.1 0.3 0; 
        0.2 0.4 0 ;
        0 0 1];
p.EGk = [0.2 0.1 0.3]';

p.H=diag([0.1 0.2 0.3]);


p.x_init = [ 1  -3   3  6   5   2  -2  -5;
             5   4   1   2  7   3  -4  -2;
             1   3  -4  3  8   4  -3  -1];
p.T = 200;
problems{3} = p;

% ===== ID 4:  =====


p = struct();
p.name = 'Exemplo 4 agentes artigo C Gain';
p.F = [sqrt(2)/2 sqrt(2)/2 0;
-sqrt(2)/2 sqrt(2)/2 0;
0 0 1];
p.Gk = [0.2 -0.4 1]';
p.Lcal = [ 2.7 -1.3 0 0 0 0 -0.7 -0.7 
-1.3 2.5 -1.2 0 0 0 0 0
 0 -1.2 2.3 -1.1 0 0 0 0
0 0 -1.1 3.1 -1.2 -0.8 0 0 
0  0 0 -1.2 2.1 -0.9 0 0
0 0 0 -0.8 -0.9 3.1 -1.4 0 
-0.7 0 0 0 0 -1.4 3.1 -1
-0.7 0 0 0 0 0 -1 1.7];

p.Ef = [0.1 0.3 0; 
        0.2 0.4 0 ;
        0 0 1];
p.EGk = [0.2 0.1 0.3]';

p.H=diag([0.1 0.2 0.3]);


p.x_init = [ 1  2   1  3   5   -3  -2  -5;
             5   4   1   2  6   3  -4  -2;
             -2   3  2  1  -2  4  -3  -1];
p.T = 120;
problems{4} = p;
% ===== ID 1:  =====
p = struct();
p.name = 'Exemplo 4 agentes artigo robusto';
p.F = [1 1;
       0 1];
p.Ef = [0.2 0.1 ;
        0.3 0.4];
p.Lcal = [ 2  -1 -1 0 0;
          -1 1 0 0 0;
           0 -1  2 -1 0;
           -1 0 -1 2 0
           0 0 0 -1 1];
p.Gk = [0.5 1]';
p.EGk = [0.2;
         0.4];
p.H=[0.2 0; 0 0.3]';
p.x_init = [ -2  4 -6 -3;
            3  -5  2 -2];
p.T = 300;
problems{5} = p;
% ===== ID 1:  =====
p = struct();
p.name = 'Exemplo 4 agentes artigo robusto';
p.F = [1 1;
       0 1];
p.Ef = [0.2 0.1 ;
        0.3 0.4];
p.Lcal = [ 2  -1 -1 0 0;
          -1 1 0 0 0;
           0 -1  2 -1 0;
           -1 0 -1 2 0
           0 0 0 -1 1];
p.Lcal = blkdiag(p.Lcal,p.Lcal*1.1,p.Lcal *0.9);
p.Gk = [0.5 1]';
p.EGk = [0.2;
         0.4];
p.H=[0.2 0; 0 0.3]';
p.x_init = [ -2  4 -6 -3;
            3  -5  2 -2];
p.T = 300;
problems{6} = p;
end
