function [Y]=adist(X,d,K)
data=X';
[m,N] = size(data);  % m is the dimensionality of the input sample points. 
% Step 0:  Neighborhood Index
M=zeros(m,m,N);
if nargin<4
    if length(K)==1
        K = repmat(K,[1,N]);
    end;
    NI = cell(1,N); 
    if m>N
        a = sum(data.*data); 
        dist2 = sqrt(repmat(a',[1 N]) + repmat(a,[N 1]) - 2*(data'*data));
        for i=1:N
            
            [dist_sort,J]=sort(dist2(:,i)); 
            Ii=J(1:K(i)); 
            NI{i} = Ii;
        end;
    else
        for i=1:N
            
            x = data(:,i); ki = K(i);
            dist2 = sum((data-repmat(x,[1 N])).^2,1);    
            [dist_sort,J] = sort(dist2);  
            Ii = J(1:ki);  
            NI{i} = Ii;
        end;
    end;
else
    K = zeros(1,N);
    for i=1:N
        K(i) = length(NI{i});
    end;
end;

%%%%%%%%%%%%%
for i=1:N
    
    Ii = NI{i}; ki = K(i);
    Xi = data(:,Ii)-repmat(mean(data(:,Ii),2),[1,ki]);
    X_mean(:,i)=mean(data(:,Ii),2);
    Xi=Xi';
    [ Y,V ]= aoge(Xi,d);
    M(:,:,i)=V*V';    
end
Iq=[];
c=zeros(1,N);
count=0;
W=zeros(m,m);
for i=1:N
 
 if c(i)==0
     Q=zeros(m,m);
    %count=count+1;
    Ii = NI{i};
    for j=Ii
    Iq=[Iq,NI{j}];
    end
    Iq=unique(Iq);
    c(Iq)=1;
    for j=Iq       
        Q=Q+M(:,:,j)*M(:,:,j)';
    end
    [U,S]=eig(Q*Q');
    
    U=U(:,1:d);    
 end
   W=W+U*U';
end

X_m=data-repmat(mean(data,2),1,N);
for i=1:N
    X_m(:,i)=X_m(:,i)/norm(X_m(:,i));
end
[U1,S]=eig(X_m*X_m'-W);
[s,Ji] = sort(-diag(S)); 
U1 = U1(:,Ji(1:d));  
Y=X*U1;
    