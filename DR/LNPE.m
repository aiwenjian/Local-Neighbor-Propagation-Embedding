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
    X = X';                                      % X.shape�� d * n
    neighborhood = neighborhood';
    neighborhood = neighborhood(2:k+1,:);        % ȡk���ڣ�ÿ�д���ͬһ����� k �����ڵ� index
    if nargout > 1
        mapping.nbhd = distance;
    end

    % Construct reconstruction weight matrix
    layerNum = 2;                                   % ���ڴ�������
    W1 = zeros(n, n, layerNum);                     % ��ά���󣬵���ά�������ֵڼ���Ȩֵ����
    X1 = zeros(d, n);                               % ���ÿ���ؽ����µĸ�ά����
    %% ��ʼ�ĵ�һ����ڣ���LLE����Ȩֵ
    for i=1:n                                       
       z = X(:,neighborhood(:,i)) - repmat(X(:,i), 1, k);		% Shift point to origin
       C = z' * z;												% Compute local covariance
       C = C + eye(k, k) * tol * trace(C);						% Regularization of covariance (if K > D)
       temp = C \ ones(k, 1);									% Solve linear system
       temp = temp / sum(temp);							% Make sure that sum is 1
       W1(neighborhood(:,i), i, 1) = temp;
       X1(:,i) = X * W1(:,i, 1);
    end
    
    %% ��2~layerNum��Ľ��ڴ���
    if layerNum >= 2
        for j = 2:layerNum                               % �� 2:layerNum ��
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
    
    % ���һ�¸�ά�����ؽ����
    disp("LLE_2layers:");
    rec_err = norm(X1-X, 'fro')
    
    %% ��άǶ��
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
