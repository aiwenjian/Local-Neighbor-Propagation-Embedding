function [mappedX, err] = Mlle(X, no_dims, k, tol)


    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 12;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end
    modified_tol=1e-12;
    % Get dimensionality and number of dimensions
    [n, D] = size(X);

    % Compute pairwise distances and find nearest neighbours (vectorized implementation)
    disp('Finding nearest neighbors...');    
    [distance, neighborhood] = find_nn(X, k + 1);
    neighborhood = neighborhood';
    neighborhood = neighborhood(2:k+1,:);
    if nargout > 1
        mapping.nbhd = distance;
    end

    % Find reconstruction weights for all points by solving the MSE problem 
    % of reconstructing a point from each neighbours. A used constraint is 
    % that the sum of the reconstruction weights for a point should be 1.
    disp('Compute reconstruction weights...');
%     if k > d 
%         tol = 1e-5;
%     else
%         tol = 0;
%     end

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
end
% calculate eta
 rho=  sum(evals(:,no_dims+1:end),2)./sum(evals(:,1:no_dims),2);
 eta=median(rho);
 
 %find s_i
   s_range=zeros(1,n);
   evals_cumsum= cumsum(evals,2);
   
  eta_range=repmat(evals_cumsum(:,end),1,nev-1)./evals_cumsum(:,1:end-1)-1;
  for i=1:n
      [~,sb]=find(eta_range(i,end:-1:1)>eta); 
      if length(sb)==0
          continue;
      else
      s_range(i)=sb(1)-1;
      end
  end
  s_range=s_range+k-nev;
  
  %calculate M
  % Now that we have the reconstruction weights matrix, we define the 
  % sparse cost matrix M = (I-W)'*(I-W).
    M = sparse(1:n, 1:n, ones(1, n), n, n, 4 * k * n);
    for i=1:n
        s_i=s_range(i);
        Vi=V(:,k-s_i+1:end,i);
        alpha_i=norm(sum(Vi,1))/sqrt(s_i);
        %compute Householder matrix 
        h=alpha_i*ones(s_i,1)-sum(Vi',2);
        norm_h= norm(h);
        if norm_h<modified_tol
            h=h*0;
        else
            h=h./norm_h;
        end
        Wi=Vi-2*(Vi*h)*h'+(1-alpha_i)*repmat(W(i,:)',1,s_i);
        
        
       w = Wi;
       j = neighborhood(:,i);
       M(j, j) = M(j, j) + w * w';
       Wi_sum=sum(Wi,2);
       M(i, j) = M(i, j) - Wi_sum';
       M(j, i) = M(j, i) - Wi_sum;
       M(i,i)=M(i,i)+s_i-1;
    end
	
	% For sparse datasets, we might end up with NaNs or Infs in M. We just set them to zero for now...
	M(isnan(M)) = 0;
	M(isinf(M)) = 0;
    % The embedding is computed from the bottom eigenvectors of this cost matrix
	disp('Compute embedding (solve eigenproblem)...');
    tol = 0;
    if strcmp(eig_impl, 'JDQR')
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [mappedX, eigenvals] = jdqr(M, no_dims + 1, tol, options);
    else
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [mappedX, eigenvals] = eigs(M, no_dims + 1, tol, options);          % only need bottom (no_dims + 1) eigenvectors
    end
    [eigenvals, ind] = sort(diag(eigenvals), 'ascend');
    if size(mappedX, 2) < no_dims + 1
		no_dims = size(mappedX, 2) - 1;
		warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
    end
    mappedX = mappedX(:,ind(2:no_dims + 1));   
    err=sum(eigenvals(2:no_dims + 1));
    % throw away zero eigenvector/value					

