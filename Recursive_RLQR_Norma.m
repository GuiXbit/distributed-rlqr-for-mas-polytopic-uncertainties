
function [P,K,L] = Recursive_RLQR_Norma(parms,mu,alfa, m, n, q, F, G, Q, R, H, Ef, Eg, P,c)

    if mu <= 0
        error('Mu must be positive')
    end 
    if (parms ~= 0) && (parms ~= 1)
        error('parms = {0,1}')
    end 
    
    In = eye(n);
    Im = eye(m);
    Il = eye(q);
    
    lambda = (1 + alfa) * norm( mu * H' * H );

    phi = mu^(-1) * In - lambda^(-1) * H * H';

    if min( eig( phi ) ) < 0
        error('ErrorTmin( eig( phi ) )ests:convertTest', 'Phi is not positive (semi)definite \n Change variables parms and/or mu.')
    end

    MSigma = parms * blkdiag( phi,lambda^(-1)*Il );
    
    Aux = [ In          -G; 
            zeros(q,n)  -Eg];

    if ( rank( Aux ) < size( Aux,1 ) ) && parms == 0
        error('Infeasible problem for parms = 0') 
    end


    P_cal = blkdiag( P,R,Q );

    A_cal = [ eye(n)      zeros(n,m); 
              zeros(m,n)  Im        ; 
              zeros(n)    zeros(n,m); 
              In         -G         ; 
              zeros(q,n) -Eg        ];

    Aux_W = blkdiag( inv(P_cal),MSigma );

    W = [ Aux_W       A_cal; 
          A_cal' zeros(n+m)];

    if rcond(W) <= eps
        error('ErrorTests:convertTest', 'Matrix W is ill-conditioned \n Change variables parms and/or mu.')
    end   

    Z = [ zeros(3*n+m+q,n) ; In ; zeros(m,n) ];

    V = [ zeros(4*n+m+q,m) ; Im ];

    U = [ zeros(n); zeros(m,n); -In; F; Ef; zeros(n+m,n) ];

    Sol = [ Z V U ]' * inv(W) * U;

    AUX = mat2cell( Sol,[n m n], n );

    L = AUX{1};

    K =-c*AUX{2};

    P = AUX{3};

    if min( eig(P) ) <= 0
        error('Riccati is not positive definite') 
    end
end

