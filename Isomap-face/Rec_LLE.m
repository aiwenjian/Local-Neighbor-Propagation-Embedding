function [Xnew] = Rec_LLE(X,K)

[D,N] = size(X);
Xnew=zeros(D,N); % 存放重建后的高维数据点
fprintf(1,'LLE running on %d points in %d dimensions\n',N,D);

% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
fprintf(1,'-->Finding %d nearest neighbours.\n',K);

X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);

% STEP2: SOLVE FOR RECONSTRUCTION WEIGHTS
fprintf(1,'-->Solving for reconstruction weights.\n');

tol=1e-5; % regularlizer in case constrained fits are ill conditioned

W = zeros(K,N);
for ii=1:N
    z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
    C = z'*z;                                        % local covariance
    C = C + eye(K,K)*tol*trace(C);                   % regularlization (K>D)
    W(:,ii) = C\ones(K,1);                           % solve Cw=1
    W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
    Xnew(:,ii) = X(:,neighborhood(:,ii)) * W(:,ii);
end

