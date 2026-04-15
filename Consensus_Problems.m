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
end
