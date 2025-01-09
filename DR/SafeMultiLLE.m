function [Y,err]=SafeMultiLLE(data,d,K,tol)

 X = data';%data(:,1:end-1)';
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
       % WW = zeros(K,K);
        for jj=1:K
            if N<400  && length(unique([neighborhood(:,neighborhood(jj,ii)),neighborhood(:,ii)]))>1.3*K                
                K1=ceil(0.5*K);
                else
                  K1=K;  
            end
            Xtemp=X(:,neighborhood(:,neighborhood(jj,ii)));
            z = Xtemp(:,1:K1)-repmat(X(:,ii),1,K1);
            C = W(jj,ii)^2*z'*z;
            C = C + eye(K1,K1)*tol*trace(C); 
            WW = C\ones(K1,1);
            WW= WW/sum(WW); 
            Ntemp=neighborhood(:,neighborhood(jj,ii));
            R(Ntemp(1:K1),jj)=WW;
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
err=sum(eigenvals(2:d + 1));

fprintf(1,'Done.\n');

