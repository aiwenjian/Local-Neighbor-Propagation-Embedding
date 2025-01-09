function [Y,T]=RHNE(data,K,d)

    X = data;%data(:,1:end-1)'; 每列是一个数据点
    [D,N] = size(X);
    X2 = sum(X.^2,1);
    tol=1e-5;
    distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
    [sorted,index] = sort(distance);
    neighborhood = index(2:(1+K),:);
    fprintf(1,'-->Solving for reconstruction weights.\n');
   % tol=1e-2;
    W = zeros(K,N);
    for ii=1:N
       z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
       C = z'*z;                                        % local covariance
            C = C + eye(K,K)*tol*trace(C); 
             % regularlization (K>D)
       W(:,ii) = C\ones(K,1);                           % solve Cw=1
       W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
    end
%     Wfull=zeros(N,N);
%     Wfull(neighborhood)=W;
    fprintf(1,'-->Solving for 2nd reconstruction weights.\n');
    
    
    %dirAML权值求解过程
    T=zeros(N,N);
    for ii=1:N
        R=zeros(N,K);
        z = X(:, neighborhood(:, neighborhood(:, ii))) - repmat(X(:,ii),1,K*K); % shift ith pt to origin
        C = z'*z;                                        % local covariance,shape(K^2, k^2)
        C = C + eye(K*K,K*K)*tol*trace(C); 
        WW = C\ones(K*K,1);                           % solve Cww=1
        WW = WW/sum(WW);                  % enforce sum(ww)=1
        reshapeWW = reshape(WW, K, K) ./ repmat(W(:, ii)', K, 1);
        sum(reshapeWW, 1)
%         reshapeWW = reshapeWW ./ repmat(sum(reshapeWW, 1), K, 1);
%         WW = reshapeWW .* repmat(W(:, ii)', K, 1);
%         WW = reshape(WW, K*K, 1);
%         sum(WW);
        for jj = 1:K
            R(neighborhood(:,neighborhood(jj,ii)),jj)=reshapeWW(:, jj);
        end
        T(:,ii)=R*W(:,ii);
    end

M2 = sparse(1:N,1:N,ones(1,N),N,N,4*K*N); 
for ii=1:N
   w = W(:,ii);
   jj = neighborhood(:,ii);
   M2(ii,jj) = M2(ii,jj) - w';
   M2(jj,ii) = M2(jj,ii) - w;
   M2(jj,jj) = M2(jj,jj) + w*w';
end;

%%Low Dimensional Embedding

% M=eye(N,N); % use a sparse matrix with storage for 4KN nonzero elements
M=(eye(N,N)-T)*(eye(N,N)-T)'+full(M2);
% CALCULATION OF EMBEDDING
options.disp = 0; options.isreal = 1; options.issym = 1; 
[Y,eigenvals] = eigs(M,d+1,0,options);
 [eigenvals, ind] = sort(diag(eigenvals), 'ascend');
Y = Y(:,ind(2:d+1))*sqrt(N); % bottom evect is [1,1,1,1...] with eval 0
Y = Y';

fprintf(1,'Done.\n');

