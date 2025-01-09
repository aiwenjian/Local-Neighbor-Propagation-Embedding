function Xnew=Rec_RHNE(X,K)

tol = 1e-4;
[D,N] = size(X);
X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);

%%%%%%%%%%%%%%%%%%  LLEȨֵ���  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'-->Solving for reconstruction weights.\n');
W = zeros(K,N);
for ii=1:N
    z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
    C = z'*z;                                        % local covariance
        C = C + eye(K,K)*tol*trace(C); 
         % regularlization (K>D)
    W(:,ii) = C\ones(K,1);                           % solve Cw=1
    W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end

%%%%%%%%%%%%%%%%%%%  dirAMLȨֵ������  %%%%%%%%%%%%%%%%%%%%
fprintf(1,'-->Solving for 2nd reconstruction weights.\n');
Xnew=zeros(D,N); % ��Ÿ�ά���ݵ�
for ii=1:N
    z = X(:, neighborhood(:, neighborhood(:, ii))) - repmat(X(:,ii),1,K*K);
    C = z'*z;                                                               % shape(K^2, K^2)
    C = C + eye(K*K,K*K)*tol*trace(C); 
    WW = C\ones(K*K,1);
    WW = WW/sum(WW);                                                        % enforce sum(ww)=1
    Xnew(:, ii) = X(:, neighborhood(:, neighborhood(:, ii))) * WW;          % ��ά�ؽ�
end