function [mappedX, err] = LNPE(X, no_dims, k, tol)

    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 6;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end

    % Get dimensionality and number of dimensions
    [n, d] = size(X);

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
    
    % 输出一下高维的总重建误差
    disp("LLE_2layers:");
    rec_err = norm(X1-X, 'fro')
    
    %% 低维嵌入
    I = eye(n);
    W = I;
    M = zeros(n, n);
    for i = 1 : layerNum
        W = W * W1(:,:,i);
        M = M + (I - W) * (I - W)';
    end

    % The embedding is computed from the bottom eigenvectors of this cost matrix
    tol = 0;
    if strcmp(eig_impl, 'JDQR')
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [mappedX, eigenvals] = jdqr(M, no_dims + 1, tol, options);
    else
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [mappedX, eigenvals] = eigs(M, no_dims + 1, tol, options);          % only need bottom (no_dims + 1) eigenvectors
    end
    [eigenvals, ind] = sort(diag(eigenvals), 'ascend');
    if size(mappedX, 2) < no_dims + 1
		no_dims = size(mappedX, 2) - 1;
		warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
    end
    mappedX = mappedX(:,ind(2:no_dims + 1));     
    err=sum(eigenvals(2:no_dims + 1));			
