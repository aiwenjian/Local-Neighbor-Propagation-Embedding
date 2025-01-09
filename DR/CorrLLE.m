function [mappedX,err]=CorrLLE(X, no_dims, k, tol)
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
W = zeros(K,N);
for ii=1:N
    z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
    C = z'*z;                                        % local covariance
    C = C + eye(K,K)*tol*trace(C);
    % regularlization (K>D)
    W(:,ii) = C\ones(K,1);                           % solve Cw=1
    W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end;
%     Wfull=zeros(N,N);
%     Wfull(neighborhood)=W;

T=zeros(N,N);

for ii=1:N
    R=zeros(N,K);
    WW = zeros(K,K);
    for jj=1:K
        z = X(:,neighborhood(:,neighborhood(jj,ii)))-repmat(X(:,ii),1,K);
        C = W(jj,ii)^2*z'*z;
        C = C + eye(K,K)*tol*trace(C);
        WW(:,jj) = C\ones(K,1);
        WW(:,jj) = WW(:,jj)/sum(WW(:,jj));
        R(neighborhood(:,neighborhood(jj,ii)),jj)=WW(:,jj);
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
M=(eye(N,N)-T)*(eye(N,N)-T)'+ full(M2);