function problems = Consensus_Problems()
% Banco unificado: cada entrada inclui dinamica + topologia + iniciais.
% Campos usados no script:
%   name, A, Lcal, Bcal, G0 ou Gk, x_init, x0_leader, T

problems = cell(1,8);

% ===== ID 1:  =====
p = struct();
p.name = 'Artigo Politopic';
A{1} = [0.75 0.75 0; -0.75 0.75 1.00; 0 0 1.05];
A{2} = [0.85 0.85 0; -0.85 0.85 0.95; 0 0 0.95];
A{3} = [0.65 0.65 0; -0.65 0.65 0.95; 0 0 0.80];
A{4} = [0.95 0.05 0; -0.05 0.95 1.00; 0 0 1.20];

B{1} = [ 1;  0;  1];
B{2} = [ 1; -1;  0];
B{3} = [-1;  0;  1];
B{4} = [ 0;  1; -1];

V = 4;

F = zeros(size(A{1}));
G = zeros(size(B{1}));
for l = 1:V
    F = F + A{l};
    G = G + B{l};
end
F = F / V;
G = G / V;

for l = 1:V
    dF{l} = (A{l} - F)*0.2;
    dG{l} = (B{l} - G)*0.2;
end

p.F = F;
p.dF = dF
p.Lcal = [ 2.7 -1.3 0 0 0 0 -0.7 -0.7 
-1.3 2.5 -1.2 0 0 0 0 0
 0 -1.2 2.3 -1.1 0 0 0 0
0 0 -1.1 3.1 -1.2 -0.8 0 0 
0  0 0 -1.2 2.1 -0.9 0 0
0 0 0 -0.8 -0.9 3.1 -1.4 0 
-0.7 0 0 0 0 -1.4 3.1 -1
-0.7 0 0 0 0 0 -1 1.7];
p.Gk = G;
p.dG = dG;
p.x_init = [ 1  2   1  3   5   -3  -2  -5;
             5   4   1   2  6   3  -4  -2;
             -2   3  2  1  -2  4  -3  -1];
p.T = 140;
problems{1} = p;
% ===== ID 2:  =====
end
