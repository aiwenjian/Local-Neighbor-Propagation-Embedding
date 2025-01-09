function Xnew = Rec_BHNE(X, K)

    tol = 1e-5;
    [D, N] = size(X);
    X2 = sum(X.^2,1);
    distance = repmat(X2, N, 1) + repmat(X2', 1, N) - 2 * X' * X;
    [sorted, index] = sort(distance);
    neighborhood = index(2:(1+K), :);

    %++++++++++++++++++++++++++++  求解LLE权值  +++++++++++++++++++++++++++
    fprintf(1,'-->Solving for reconstruction weights.\n');
    W = zeros(K,N);
    for ii=1:N
        z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K);     % shift ith pt to origin
        C = z'*z;                                            % local covariance
        C = C + eye(K,K)*tol*trace(C);                       % regularlization (K>D)     
        W(:,ii) = C\ones(K,1);                               % solve Cw=1
        W(:,ii) = W(:,ii)/sum(W(:,ii));                      % enforce sum(w)=1
    end

    % ++++++++++++++++++++++++++ BHNE权值求解过程  ++++++++++++++++++++++
    fprintf(1,'-->Solving for 2nd reconstruction weights.\n');
    tol=1e-6;
    EPOCH = 2;  % 迭代轮次
    Xnew=zeros(D,N);
    for ii=1:N
        %第一轮得到初始值，从第二轮迭代开始进行更新
        WW = zeros(K*K, 1);
        for epoch = 1:EPOCH
            for j = 1:K                                                              %j是Xi的第j个近邻点
                if epoch == 1                                                        %第一轮基于LLE
                    W_i = W(:,ii) ; W_i(j)=0;
                    recX_i = X(:, neighborhood(:, ii)) * W_i;
                    X_ij = X(:,ii) - recX_i;                                           %计算要更新的点值
                else
                    recX_i = X(:, neighborhood(:, neighborhood(:, ii))) * WW;            % recX_i是利用直接权值求解算出来的重建X_i
                    ori_WW = WW((j-1)*K+1:j*K);                                          %该轮需要更新的WW权值
                    objX_ij = X(:, neighborhood(:, neighborhood(j, ii))) * ori_WW;
                    X_ij = X(:,ii) - (recX_i - objX_ij);
                end  
                X_ij_ = W(j, ii) * X(:, neighborhood(:, neighborhood(j, ii))) - repmat(X_ij, 1, K);
                CC = X_ij_' * X_ij_;
                CC = CC + eye(K, K) * tol * trace(CC);
                W_ij = CC \ ones(K, 1);
                W_ij = W_ij / sum(W_ij);                                                %更新完的第二层权值
                WW((j-1)*K+1:j*K) = W_ij * W(j, ii);
            end
            recX_i = X(:, neighborhood(:, neighborhood(:, ii))) * WW;
        end
        Xnew(:, ii) = X(:, neighborhood(:, neighborhood(:, ii))) * WW;
    end