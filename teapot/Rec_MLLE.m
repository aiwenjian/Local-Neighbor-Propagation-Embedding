function Xnew = Rec_MLLE(X, k)

[D,n] = size(X);
X = X';
Xnew=zeros(D,n); % 存放重建后的高维数据点

% Compute pairwise distances and find nearest neighbours (vectorized implementation)
disp('Finding nearest neighbors...');    
[distance, neighborhood] = find_nn(X, k + 1);
neighborhood = neighborhood';
neighborhood = neighborhood(2:k+1,:);
if nargout > 1
    mapping.nbhd = distance;
end

disp('Compute reconstruction weights...');
% Construct reconstruction weight matrix
V = zeros(k,k,n);
nev = min(D,k);
evals = zeros(n,nev);

if k>D
    for i=1:n
        z = X(neighborhood(:,i),:) - repmat(X(i,:), k, 1);
        [V(:,:,i),eigvalues,~] = svd(z);
        evals(i,:)=diag(eigvalues);
    end
    evals=evals.^2;
else
    for i=1:n
        z = X(neighborhood(:,i),:) - repmat(X(i,:), k, 1); 
        C = z * z';	
        [vi,eivalue]=eig(C);
        Tempeig=diag(eivalue);
        evals(i,:)=Tempeig(end:-1:1);
        V(:,:,i)=vi(:,end:-1:1);
    end
end

reg=(1e-3)*sum(evals,2);

tmp1=permute(V,[2,1,3]);
tmp1=sum(tmp1,2);
tmp=[reshape(tmp1,k,n)]';

tmp(:,1:nev)=tmp(:,1:nev)./(evals+repmat(reg,1,nev));
tmp(:,nev+1:end)=tmp(:,nev+1:end)./repmat(reg,1,k-nev);

W = zeros(n,k);
for i=1:n					
    W(i,:) = tmp(i,:)*[V(:,:,i)]';
    W(i,:) = W(i,:) / sum(W(i,:));	
    Xnew(:, i) = X(neighborhood(:, i), :)' * W(i, :)'; % 高维重建
end
