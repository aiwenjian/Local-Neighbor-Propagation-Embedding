function [Xnew] = Rec_LNPE(X, k)

    % Get dimensionality and number of dimensions
    X = X';
    [n, d] = size(X);
    tol = 1e-5;

    % Compute pairwise distances and find nearest neighbours (vectorized implementation)  
    [distance, neighborhood] = find_nn(X, k + 1);
    X = X';                                      % X.shape： d * n
    neighborhood = neighborhood';
    neighborhood = neighborhood(2:k+1,:);        % 取k近邻，每列代表同一个点的 k 个近邻的 index
    if nargout > 1
        mapping.nbhd = distance;
    end

    % Construct reconstruction weight matrix
    layerNum = 2;                                   % 近邻传播次数
    W1 = zeros(n, n, layerNum);                     % 三维矩阵，第三维用来区分第几层权值矩阵
    X1 = zeros(d, n);                               % 存放每层重建后新的高维数据
    %% 初始的第一层近邻，用LLE计算权值
    for i=1:n                                       
       z = X(:,neighborhood(:,i)) - repmat(X(:,i), 1, k);		% Shift point to origin
       C = z' * z;												% Compute local covariance
       C = C + eye(k, k) * tol * trace(C);						% Regularization of covariance (if K > D)
       temp = C \ ones(k, 1);									% Solve linear system
       temp = temp / sum(temp);							% Make sure that sum is 1
       W1(neighborhood(:,i), i, 1) = temp;
       X1(:,i) = X * W1(:,i, 1);
    end
    
    %% 第2~layerNum层的近邻传播
    if layerNum >= 2
        for j = 2:layerNum                               % 第 2:layerNum 次
            for i=1:n                       
                z = X1(:,neighborhood(:,i)) - repmat(X(:,i), 1, k);
                C = z' * z;
                C = C + eye(k, k) * tol * trace(C);
                temp = C \ ones(k, 1);
                temp = temp / sum(temp);	
                W1(neighborhood(:,i), i, j) = temp;
                X1(:,i) = X1 * W1(:,i, j);
            end
        end
    end
    
    Xnew = X1;
   