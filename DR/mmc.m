function [Y,eigvector]= mmc(Z,gnd,Dim)

% sig=max(max(Z));%(max(max(Z))-min(min(Z)))/
% Z=Z/sig;
% ====== Initialization
[nSmp,D] = size(Z);
if length(gnd) ~= nSmp
    error('gnd and data mismatch!');
end

classLabel = unique(gnd);
nClass = length(classLabel);
Xa=[];
T_time=cputime;
for i = 1:nClass
    index = find(gnd==classLabel(i));
    classMean(i,:) = mean(Z(index,:),1);
    X=Z(index,:);
    if length(index)==1
        Xc(1,:)=X(1,:)/norm(X(1,:));
    else
        for j=1:length(index)
            Xc(j,:)=X(j,:)-classMean(i,:);
        end        
    end    
   Xa=[Xa;Xc];
end


Xb=zeros(length(nClass),D);

globalMean=mean(classMean,1);
for i = 1:nClass 
     index = find(gnd==classLabel(i));
    Xb(i,:)=sqrt(length(index))*(classMean(i,:)-globalMean);
end

Xg=Xb'*Xb;
Xw=-Xa'*Xa;


    M=Xg+Xw;
     M=(M+M')/2;
    [U,SS]=eig(M);%jdqr(M,Dim);
    eigvalue=diag(SS);
    
    [S,ind]=sort(eigvalue,'descend');
    eigvector=U(:,ind(1:Dim));
    Y=Z*eigvector;
    mmctime=cputime-T_time;


