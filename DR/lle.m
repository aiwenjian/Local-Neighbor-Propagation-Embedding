function [mappedX, err] = lle(X, no_dims, k, tol)
    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 12;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end

    % Get dimensionality and number of dimensions
    [n, d] = size(X);

    % Compute pairwise distances and find nearest neighbours (vectorized implementation)
    disp('Finding nearest neighbors...');    
    [distance, neighborhood] = find_nn(X, k + 1);
    X = X';
    neighborhood = neighborhood';
    neighborhood = neighborhood(2:k+1,:);
    if nargout > 1
        mapping.nbhd = distance;
    end

    % Find reconstruction weights for all points by solving the MSE problem 
    % of reconstructing a point from each neighbours. A used constraint is 
    % that the sum of the reconstruction weights for a point should be 1.
    disp('Compute reconstruction weights...');


    % Construct reconstruction weight matrix
    W = zeros(k, n);
    for i=1:n
       z = X(:,neighborhood(:,i)) - repmat(X(:,i), 1, k);		% Shift point to origin
       C = z' * z;												% Compute local covariance
       C = C + eye(k, k) * tol * trace(C);						% Regularization of covariance (if K > D)
       W(:,i) = C \ ones(k, 1);									% Solve linear system
       W(:,i) = W(:,i) / sum(W(:,i));							% Make sure that sum is 1
    end

    % Now that we have the reconstruction weights matrix, we define the 
    % sparse cost matrix M = (I-W)'*(I-W).
    M = sparse(1:n, 1:n, ones(1, n), n, n, 4 * k * n);
    for i=1:n
       w = W(:,i);
       j = neighborhood(:,i);
       M(i, j) = M(i, j) - w';
       M(j, i) = M(j, i) - w;
       M(j, j) = M(j, j) + w * w';
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

