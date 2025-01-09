function [Y,eigvector]= LDAE(Z,gnd,Dim)

% sig=max(max(Z));%(max(max(Z))-min(min(Z)))/
% Z=Z/sig;
% ====== Initialization
[nSmp,D] = size(Z);

if length(gnd) ~= nSmp
    error('gnd and data mismatch!');
end

classLabel = unique(gnd);
nClass = length(classLabel);
Xw=zeros(D,D);
T_time=cputime;
for i = 1:nClass
    index = find(gnd==classLabel(i));
    classMean(i,:) = mean(Z(index,:),1);
    X=Z(index,:);
    Xc=zeros(length(index),D);
    if length(index)==1
        Xc(1,:)=X(1,:)/norm(X(1,:));
    else
        for j=1:length(index)
            Xc(j,:)=X(j,:)-classMean(i,:);
            Xc(j,:)=Xc(j,:)/norm(Xc(j,:));
        end
    end
        Xw=Xw-Xc'*Xc;

end
Xb=zeros(length(nClass),D);
globalMean=mean(classMean,1);
for i = 1:nClass
    index = find(gnd==classLabel(i));
    Xb(i,:)=classMean(i,:)-globalMean;
    
    a(i)=length(index);
    Xb(i,:)=sqrt(a(i))*Xb(i,:)/norm(Xb(i,:));
    
end

    Xg=Xb'*Xb;
    for i=1:D
        Xw(i,i)=nSmp/Dim+Xw(i,i);
    end
    M=Xg+Xw;%max
    M=(M+M')/2;
    options.LSolver = 'bicgstab';
    [U,SS]=jdqr(M,Dim);%eig(M);%lanczos(M,Dim);
    eigvalue=diag(SS);
    [S,ind]=sort(eigvalue,'descend');
    eigvector=U(:,ind(1:Dim));
    Y=Z*eigvector;
atime=cputime-T_time;
%save ALDE  M;


