function [Xnew]=Rec_IHNE(data,K)

X = data;
[D,N] = size(X);

X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);
fprintf(1,'-->Solving for reconstruction weights.\n');
tol=1e-4;
W = zeros(K,N);
for ii=1:N
    z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
    C = z'*z;                                        % local covariance
    C = C + eye(K,K)*tol*trace(C); 
    W(:,ii) = C\ones(K,1);                           % solve Cw=1
    W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end

Xnew=zeros(D,N); % 存放重建后的高维数据点

T=zeros(N,N);
for ii=1:N
    R=zeros(N,K);
    WW = zeros(K,K);
    TempXi=[];
    for jj=1:K
        z = X(:,neighborhood(:,neighborhood(jj,ii)))-repmat(X(:,ii),1,K);
        C = W(jj,ii)^2*z'*z;
        C = C + eye(K,K)*tol*trace(C); 
        WW(:,jj) = C\ones(K,1);
        WW(:,jj) = WW(:,jj)/sum(WW(:,jj)); 
        R(neighborhood(:,neighborhood(jj,ii)),jj)=WW(:,jj);
        TempXi=[TempXi,sum((ones(D,1)*[WW(:,jj)]').*X(:,neighborhood(:,neighborhood(jj,ii))),2)];
    end
    Xnew(:,ii)=sum((ones(D,1)*[W(:,ii)]').*TempXi,2); % 高维重建
end