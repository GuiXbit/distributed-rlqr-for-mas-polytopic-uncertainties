
function [P,K,L] = Recursive_LQR_c(parms, mu, m, n, F, G, Q, R, P,c)
    if mu <= 0
     error('Mu must be positive')
    end
    
    if (parms ~= 0) && (parms ~= 1)
      error('parms = {0,1}')
    end
    
    In = eye(n);
    Im = eye(m);
    
    Aux = [ In -G ];
    if ( rank( Aux ) < size( Aux,1 ) ) && parms ==0
      error('Infeasible problem for parms = 0') 
    end
    
    MSigma = parms * mu^(-1) * In;
                  
    P_cal = blkdiag(P,R,Q);
    
    A_cal = [ eye(n)      zeros(n,m);
              zeros(m,n)  Im        ; 
              zeros(n)    zeros(n,m); 
              In         -G         ];
    
    Aux_W = blkdiag( inv(P_cal),MSigma );
    
    W = [ Aux_W  A_cal; A_cal' zeros(n+m) ];
    
    if rcond(W) <= eps
       error('ErrorTests:convertTest','Matrix W is ill-conditioned \n Change variables parms and/or mu.')
    end
    
    Z = [ zeros(3*n+m,n) ; In ; zeros(m,n) ];
    
    V = [ zeros(4*n+m,m) ; Im ];
    
    U = [ zeros(n); zeros(m,n); -In; F; zeros(n+m,n) ]; 
    
    Sol = [ Z V U ]' * inv(W) * U;
    
    Aux = mat2cell( Sol,[n m n], n );
    
    L = Aux{1};
    
    K = -c*Aux{2};
    
    P = Aux{3};
    
    if min( eig( P ) ) <= 0
        error('Riccati is not positive definite') 
    end
end 
