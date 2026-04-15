function [P,K,L] = Recursive_Lqr_Consensus(F,G,Q,R,P,delta)

    % Verifica se as matrizes F e G são compatíveis
    if size(F, 1) ~= size(G, 1)
        error('As matrizes F e G devem ter o mesmo número de linhas.');
    end
    
    % Verifica se a matriz P é quadrada e compatível com F
    if size(P, 1) ~= size(F, 1) || size(P, 2) ~= size(F, 1)
        error('A matriz P deve ser quadrada e compatível com F.');
    end
    
    
    % Calcula o ganho K usando a fórmula do LQR
    K = inv(G' * P * G + R) * (G' * P * F);
    L=0;
    % Atualiza a matriz P usando a equação de Riccati
    P = F'*P*F-(1-delta^2)*F'*P*G*inv(G'*P*G + R)*G'*P*F +Q;
end
