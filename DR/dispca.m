function [Y,Q]=dispca(data,d,K)
data=data';
[m,N] = size(data);  % m is the dimensionality of the input sample points. 
% Step 0:  Neighborhood Index
if nargin<4
    if length(K)==1
        K = repmat(K,[1,N]);
    end;
    NI = cell(1,N); 
    if m>N
        a = sum(data.*data); 
        dist2 = sqrt(repmat(a',[1 N]) + repmat(a,[N 1]) - 2*(data'*data));
        for i=1:N
            % Determine ki nearest neighbors of x_j
            [dist_sort,J]=sort(dist2(:,i)); 
            Ii=J(1:K(i)); 
            NI{i} = Ii;
        end;
    else
        for i=1:N
            % Determine ki nearest neighbors of x_j
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
% Step 1:  local information
U=zeros(m,m);
for i=1:N
    % Compute the d largest right singular eigenvectors of the centered matrix
    Ii = NI{i}; ki = K(i);
    Xi = data(:,Ii)-repmat(mean(data(:,Ii),2),[1,ki]);  
    for j=1:ki
        Xii(:,j)=Xi(:,j)/norm(Xi(:,j));
    end
    [Vi,Si]=effPCA(Xii,d);
    [s,Ji] = sort(-diag(Si)); 
    Vi = Vi(:,Ji(1:d));  
    U=U+Vi*Vi';%*Xi*Xi';
    
end
W=data*(eye(N)-1/N*ones(N,1)*ones(1,N))*data'*U;
W=(W+W')/2;
[Q,S]=eig(W);
%[Q,S,Q1]=svd(data*(eye(N)-1/N*ones(N,1)*ones(1,N))*data'*U);
%Q=Q(:,end-d+1:end);
diag(S)
Q=Q(:,1:d);
Y=data'*Q;