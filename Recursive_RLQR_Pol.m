function [P_rlqr,K,L,Nf_hat,Ng_hat] = Recursive_RLQR_Pol(parms,mu,beta,n, m, v, F, G_rlqr, Nf, Ng, Q, R, P_rlqr)

    if mu <= 0
        error('Mu must be positive')
    end    
    if (parms ~= 0) && (parms ~= 1)
        error('parms = {0,1}')
    end 
    
    In = eye(n);
    Im = eye(m);
    Iv = eye(v*n);
    
    lambda = (beta * mu)^(-1) * Iv; 
    
    AUX1 = [];
    AUX2 = [];
    AUX3 = [];
    AUX4 = [];
    AUX5 = [];
    
    for vert = 1:v
        AUX1 = [AUX1; F];
        AUX2 = [AUX2; G_rlqr];
        AUX3 = [AUX3; Nf(:,:,vert)];
        AUX4 = [AUX4; Ng(:,:,vert)];
        AUX5 = [AUX5; eye(n)];
    end
     
    F_hat  = AUX1;
    G_hat  = AUX2;
    Nf_hat = AUX3;
    Ng_hat = AUX4;
    I_hat  = AUX5;
    
    Aux = [ I_hat -G_hat ; zeros(n*v,n) -Ng_hat ];
    if ( rank( Aux ) < size( Aux,1 ) ) && parms == 0
        error('Infeasible problem for parms = 0') 
    end
    
    phi =  mu^(-1) * Iv - lambda;    
    
    if min( eig( phi ) ) < 0
        error('ErrorTmin( eig( phi ) )ests:convertTest', 'Phi is not positive (semi)definite \n Change variables parms and/or mu.')
    end
    
    MSigma = parms * blkdiag( phi,lambda );
    norm(MSigma);
    
    P_cal = blkdiag( P_rlqr,R,Q );
                  
    A_cal = [ eye(n)        zeros(n,m); ...
            zeros(m,n)    Im        ; ...
            zeros(n)      zeros(n,m); ...
            I_hat         -G_hat     ; ...
            zeros(v*n,n)  -Ng_hat        ];
          
    
    Aux_W = blkdiag( inv(P_cal),MSigma );
            
    W = [ Aux_W  A_cal; A_cal' zeros(n+m) ];
    
    
    if rcond(W) <= eps
        error('ErrorTests:convertTest', 'Matrix W is ill-conditioned \n Change variables parms and/or mu.')
    end     
    
    Z = [ zeros(2*n+m+2*v*n,n) ; In ; zeros(m,n) ];
    
    V = [ zeros(3*n+m+2*v*n,m) ; Im ];
    
    U = [ zeros(n); zeros(m,n); -In; F_hat; Nf_hat; zeros(n+m,n) ];
    
    
    Sol = [ Z V U ]' * inv(W) * U;
    
    AUX = mat2cell( Sol,[n m n],[n] );
    
    L = AUX{1};
    
    K = AUX{2};
    
    P_rlqr = AUX{3};
    
    if min( eig( P_rlqr ) ) <= 0
        error('Riccati is not positive definite') 
    end
end
