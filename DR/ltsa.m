function [mappedX,err] = ltsa(X, no_dims, k, eig_impl)
%LTSA Runs the local tangent space alignment algorithm
%
%   mappedX = ltsa(X, no_dims, k, eig_impl)
%
% The function runs the local tangent space alignment algorithm on dataset
% X, reducing the data to dimensionality d. The number of neighbors is
% specified by k.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.2b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want. However, it is appreciated if you maintain the name of the original
% author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 12;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end
 
    % Compute neighborhood indices
    disp('Find nearest neighbors...');
    n = size(X, 1);
    [D, ni] = find_nn(X, k); 

    % Compute local information matrix for all datapoints
    disp('Compute local information matrices for all datapoints...');
    Bi = {}; 
    for i=1:n
        % Compute correlation matrix W
        Ii = ni(i,:);
        Xi = X(Ii,:) - repmat(mean(X(Ii,:), 1), [k 1]);
        W = Xi * Xi'; 
        W = (W + W') / 2;
        % Compute local information by computing d largest eigenvectors of W
        [Vi, Si] = schur(W);
        [s, Ji] = sort(-diag(Si));
		if length(Ji) < no_dims
			no_dims = length(Ji);
			warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
		end
        Vi = Vi(:,Ji(1:no_dims));  
        % Store eigenvectors in G (Vi is the space with the maximum variance, i.e. a good approximation of the tangent space at point Xi)
		% The constant 1/sqrt(k) serves as a centering matrix
		Gi = double([repmat(1 / sqrt(k), [k 1]) Vi]);
		% Compute Bi = I - Gi * Gi'
		Bi{i} = eye(k) - Gi * Gi';  
    end
    
    % Construct sparse matrix B (= alignment matrix)
    disp('Construct alignment matrix...');
    B = speye(n);
    for i=1:n
        Ii = ni(i,:);
        B(Ii, Ii) = B(Ii, Ii) + Bi{i};							% sum Bi over all points
		B(i, i) = B(i, i) - 1;
    end
	B = (B + B') / 2;											% make sure B is symmetric
	
	% For sparse datasets, we might end up with NaNs in M. We just set them to zero for now...
	B(isnan(B)) = 0;
	B(isinf(B)) = 0;
    
    % Perform eigenanalysis of matrix B
	disp('Perform eigenanalysis...');
    tol = 0;
	if  condest(B)==inf%strcmp(eig_impl, 'JDQR')
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [mappedX, D] = jdqr(B, no_dims + 2, tol, options);      % only need bottom (no_dims + 1) eigenvectors
    else
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [mappedX, D] = eigs(B, no_dims + 2, tol, options);      % only need bottom (no_dims + 1) eigenvectors
    end

    % Sort eigenvalues and eigenvectors
    [evals, ind] = sort(diag(D), 'ascend');
    mappedX = mappedX(:,ind);

    % Final embedding coordinates
	if size(mappedX, 2) < no_dims + 1, no_dims = size(mappedX, 2) - 1; end
    mappedX = mappedX(:,2:no_dims + 1);
    err=sum(evals(2:no_dims + 1));