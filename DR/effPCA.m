function [V,Si]=effPCA(X1,d)
[m,n]=size(X1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%LOSA中
 X=X1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%真正的Eff
% 
% X=X1*(eye(n,n)-1/n*ones(n,1)*ones(1,n));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W=X'*X; W=(W+W')/2;
[U,Si]=eig(W);
lambda=diag(Si);
V=zeros(m,size(U,2));
for j=1:size(U,2)
    V(:,j)=(1/lambda(j))*X*U(:,j);
end
