function [mappedX,err]=CorrLNPE(X, no_dims, k, tol)
eig_impl= 'matlab';
n = size(X, 1);
% index=randperm(n);
% X=X(index,:);
nCor=ceil(0.4*n);
CorrX=X(1:nCor,:);
Xs=X(nCor+1:n,:);
X1=Xs(1:floor((n-nCor)/2),:);
X2=Xs(floor((n-nCor)/2)+1:end,:);
%%
M1=mllew([CorrX;X1],no_dims,k, tol);
M2=mllew([CorrX;X2],no_dims,k, tol);
M1cc=M1(1:nCor,1:nCor);
M1cs=M1(1:nCor,nCor+1:end);
M1sc=M1(nCor+1:end,1:nCor);
M1ss=M1(nCor+1:end,nCor+1:end);

%%
M2cc=M2(1:nCor,1:nCor);
M2cs=M2(1:nCor,nCor+1:end);
M2sc=M2(nCor+1:end,1:nCor);
M2ss=M2(nCor+1:end,nCor+1:end);

MM=[M1cc+M2cc M1cs M2cs;M1sc M1ss zeros(size(M1ss,1),size(M2cs,2));
    M2sc zeros(size(M2sc,1),size(M1ss,2)) M2ss];
disp('Compute embedding (solve eigenproblem)...');
tol = 0;
if strcmp(eig_impl, 'JDQR')
    options.Disp = 0;
    options.LSolver = 'bicgstab';
    [mappedX, eigenvals] = jdqr(MM, no_dims + 1, tol, options);
else
    options.disp = 0;
    options.isreal = 1;
    options.issym = 1;
    [mappedX, eigenvals] = eigs(MM, no_dims + 1, tol, options);          % only need bottom (no_dims + 1) eigenvectors
end
[eigenvals, ind] = sort(diag(eigenvals), 'ascend');
if size(mappedX, 2) < no_dims + 1
    no_dims = size(mappedX, 2) - 1;
    warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
end
mappedX = mappedX(:,ind(2:no_dims + 1));
err=sum(eigenvals(2:no_dims + 1));

% function [M]=mllew(X, no_dims, k, tol)
% disp('Finding nearest neighbors...');
% [distance, neighborhood] = find_nn(X, k + 1);
% X = X';
% n=size(X,2);
% neighborhood = neighborhood';
% neighborhood = neighborhood(2:k+1,:);
% if nargout > 1
%     mapping.nbhd = distance;
% end
% disp('Compute reconstruction weights...');
% % if k >  no_dims
% %     tol = 1e-4;
% % else
% %     tol = 0;
% % end
% 
% % Construct reconstruction weight matrix
% W = zeros(k, n);
% for i=1:n
%     z = X(:,neighborhood(:,i)) - repmat(X(:,i), 1, k);		% Shift point to origin
%     C = z' * z;												% Compute local covariance
%     C = C + eye(k, k) * tol * trace(C);						% Regularization of covariance (if K > D)
%     W(:,i) = C \ ones(k, 1);									% Solve linear system
%     W(:,i) = W(:,i) / sum(W(:,i));							% Make sure that sum is 1
% end
% 
% % Now that we have the reconstruction weights matrix, we define the
% % sparse cost matrix M = (I-W)'*(I-W).
% M = sparse(1:n, 1:n, ones(1, n), n, n, 4 * k * n);
% for i=1:n
%     w = W(:,i);
%     j = neighborhood(:,i);
%     M(i, j) = M(i, j) - w';
%     M(j, i) = M(j, i) - w;
%     M(j, j) = M(j, j) + w * w';
% end
% 
% % For sparse datasets, we might end up with NaNs or Infs in M. We just set them to zero for now...
% M(isnan(M)) = 0;
% M(isinf(M)) = 0;

function  [M]=mllew(X, no_dims, K, tol)
tol=1e-4;
X=X';
[D,N] = size(X);
X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);
fprintf(1,'-->Solving for reconstruction weights.\n');
%tol=1e-2;

layerNum = 3;  
W1 = zeros(N, N, layerNum);                     % 三维矩阵，第三维用来区分第几层权值矩阵
X1 = zeros(D, N);

for ii=1:N   
    z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
    C = z'*z;                                        % local covariance
    C = C + eye(K,K)*tol*trace(C);
    % regularlization (K>D)
    temp = C \ ones(K, 1);
    temp = temp / sum(temp);
    W1(neighborhood(:,ii), ii, 1) = temp;                           % solve Cw=1
    X1(:,ii) = X * W1(:,ii, 1);                  % enforce sum(w)=1
end;

 %% 第2~layerNum层的近邻传播
    if layerNum >= 2
        for j = 2:layerNum                               % 第 2:layerNum 次
            for i=1:N                       
                z = X1(:,neighborhood(:,i)) - repmat(X(:,i), 1, K);
                C = z' * z;
                C = C + eye(K, K) * tol * trace(C);
                temp = C \ ones(K, 1);
                temp = temp / sum(temp);	
                W1(neighborhood(:,i), i, j) = temp;
                X1(:,i) = X1 * W1(:,i, j);  
            end
        end
    end
    
%     Wfull=zeros(N,N);
%     Wfull(neighborhood)=W;

    I = eye(N);
    W = I;
    M2 = zeros(N, N);
    for i = 1 : layerNum
        W = W * W1(:,:,i);
        M2 = M2 + (I - W) * (I - W)';
    end
%%Low Dimensional Embedding

% M=eye(N,N); % use a sparse matrix with storage for 4KN nonzero elements
M = M2;
