function XSim = simulateOptimalReach(A, A_hat, B, B_hat, C, C0, H, H_hat, D, D0, E0, X1, S1, L, K, N, FBtype)

szX = size(A,1);
szY = size(H,1);
szC = size(C,3);
szD = size(D,3);

% square root of S1
[u,s,v] = svd(S1);
sqrtS = u*diag(sqrt(diag(s)))*v';
XSim = zeros(szX,N);
Xhat = zeros(szX,N);
Xstar = zeros(szX,N);
XSim(:,1) = X1 + sqrtS*randn(szX,1);
Xhat(:,1) = X1;



for k=1:N-1
    % Update control 
    U = -L(:,:,k)*Xhat(:,k);
    
    % Compute multiplicative noise term C
    mntC = 0;
    for i=1:szC
        mntC = mntC + (C(:,:,i)*U).*randn([szX 1]); 
    end
    

    % Compute next state
    XSim(:,k+1)   = A*XSim(:,k) + B*U + C0*randn(szX,1) + mntC;

    % Predict next state
    CorrDis = A_hat*Xhat(:,k) + B_hat*U + E0*randn(szX,1);
    Xstar(:,k+1) = CorrDis;

    

    switch FBtype
        case 'Todorov2005'
            % Compute current multiplicative noise term D
            mntD = 0;
            for i=1:szD
                mntD = mntD + (D(:,:,i)*XSim(:,k)).*randn([szY 1]); 
            end

            % Compute current observation
            y = H*XSim(:,k) + D0*randn(szY,1) + mntD;

            % Estimate next state
            Xhat(:,k+1) = CorrDis + K(:,:,k)*(y-H_hat*Xhat(:,k));


        case 'Crevecoeur2011'
            % Compute next multiplicative noise term D
            mntD = 0;
            for i=1:szD
                mntD = mntD + (D(:,:,i)*XSim(:,k+1)).*randn([szY 1]); 
            end

            % Compute next observation
            y = H*XSim(:,k+1) + D0*randn(szY,1) + mntD;

            % Estimate next state
            Xhat(:,k+1) = CorrDis + K(:,:,k)*(y-H_hat*CorrDis);
    end


end