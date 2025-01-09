function [mappedA] = compute_mapping(A, type, no_dims, varargin)
%COMPUTE_MAPPING Performs dimensionality reduction on a dataset
%
%   mappedA = compute_mapping(A, type)
%   mappedA = compute_mapping(A, type, no_dims)
%
% Performs a technique for dimensionality reduction on the data specified
% in A, reducing data with a lower dimensionality in mappedA.
% The data on which dimensionality reduction is performed is given in A
% (rows correspond to observations, columns to dimensions). A may also be a
% (labeled or unlabeled) PRTools dataset.
% The type of dimensionality reduction used is specified by type. Possible
% values are 'PCA', 'SPCA', 'LDA', 'ICA', 'MDS', 'Isomap', 'LandmarkIsomap',
% 'LLE', 'Laplacian', 'HessianLLE', 'LTSA', 'DiffusionMaps', 'KernelPCA',
% 'GDA', 'SNE', 'SPE', 'AutoEncoder', and 'AutoEncoderEA'.
% The function returns the low-dimensional representation of the data in the
% matrix mappedA. If A was a PRTools dataset, then mappedA is a PRTools
% dataset as well. For some techniques, information on the mapping is
% returned in the struct mapping.
% The variable no_dims specifies the number of dimensions in the embedded
% space (default = 2). For 'LDA' and 'GDA', the labels of the instances
% should be specified in the first column of A.
%
%   mappedA = compute_mapping(A, type, no_dims, parameters)
%   mappedA = compute_mapping(A, type, no_dims, parameters, eig_impl)
%
% Free parameters of the techniques can be defined as well (on the place of
% the dots). These parameters differ per technique, and are listed below.
% For techniques that perform spectral analysis of a sparse matrix, one can
% also specify in eig_impl the eigenanalysis implementation that is used.
% Possible values are 'Matlab' and 'JDQR' (default = 'Matlab'). We advice
% to use the 'Matlab' for datasets of with 10,000 or less datapoints;
% for larger problems the 'JDQR' might prove to be more fruitful.
% The free parameters for the techniques are listed below (the parameters
% should be provided in this order):
%   LOSA  这个是当前要研究的主要算法，这个算法的特征分解矩阵
%   的奇异可能程度比较大，主要由cond（B）来判断。
%   LLOSA 通过变换得到的线性的LOSA，用于识别问题
%   PCA:            - none
%   LDA:            - none
%   ICA:            - none
%   MDS:            - none
%   Isomap:         - <int> k -> default = 12
%   LandmarkIsomap: - <int> k -> default = 12
%                   - <double> percentage -> default = 0.2
%   LLE:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   LLC:            - <int> k -> default = 12
%                   - <int> no_analyzers -> default = 20
%                   - <int> max_iterations -> default = 200
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   Laplacian:      - <int> k -> default = 12
%                   - <double> sigma -> default = 1.0
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   HessianLLE:     - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   LTSA:           - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   MVU:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   CCA:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   FastMVU:        - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   DiffusionMaps:  - <double> sigma -> default = 1.0
%   KernelPCA:      - <char[]> kernel -> {'linear', ['gauss'], 'poly', 'subsets'}
%                   - kernel parameters: type HELP GRAM for info
%   GDA:            - <char[]> kernel -> {'linear', ['gauss'], 'poly', 'subsets'}
%                   - kernel parameters: type HELP GRAM for info
%   SNE:            - <double> sigma -> default = 1.0
%   LPP:            - <int> k -> default = 12
%   NPE:            - <int> k -> default = 12
%   LLTSA:          - <int> k -> default = 12
%   SPE:            - <char[]> type -> {['Global'], 'Local'}
%                   - if 'Local': <int> k -> default = 12
%   AutoEncoder:    - none
%   AutoEncoderEA:  - none
%
% In the parameter list above, {.., ..} indicates a list of options, and []
% indicates the default setting. The variable k indicates the number of
% nearest neighbors in a neighborhood graph. The variable sigma indicates
% the variance of a Gaussian kernel.

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.2b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want. However, it is appreciated if you maintain the name of the original
% author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

% Check inputs
if nargin < 2
    error('Function requires at least two inputs.');
end
if ~exist('no_dims', 'var')
    no_dims = 2;
end
if ~isempty(varargin) && strcmp(varargin{length(varargin)}, 'JDQR')
    eig_impl = 'JDQR';
    varargin(length(varargin)) = [];
elseif ~isempty(varargin) && strcmp(varargin{length(varargin)}, 'Matlab')
    eig_impl = 'Matlab';
    varargin(length(varargin)) = [];
else
    eig_impl = 'Matlab';
end
mapping = struct;

% Handle PRTools dataset
if strcmp(class(A), 'dataset')
    prtools = 1;
    AA = A;
    if ~strcmp(type, {'LDA', 'FDA', 'GDA', 'KernelLDA', 'KernelFDA'})
        A = A.data;
    else
        A = [double(A.labels) A.data];
    end
else
    prtools = 0;
end

% Make sure we are working with doubles
A = double(A);

% Check whether value of no_dims is correct
if no_dims < 1 || no_dims > size(A, 2) || round(no_dims) ~= no_dims
    error('Value of no_dims should be a positive integer smaller than the original data dimensionality.');
end

% Switch case
switch type     
    case {'HLLE', 'HessianLLE'}
        if isempty(varargin), mappedA = hlle(A, no_dims,varargin{1},varargin{2});
        else mappedA = hlle(A, no_dims, varargin{1}, eig_impl); end
        mapping.name = 'HLLE';
        
    case 'LLE'
        if length(varargin)>1, [mappedA,err ]= lle(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err]= lle(A, no_dims, varargin{1}); end
        mapping.name = 'LLE';
        
    case 'CorrLLEO'
        if length(varargin)>1, [mappedA,err ]= CorrLLEO(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err]= CorrLLEO(A, no_dims, varargin{1}); end
        mapping.name = 'CorrLLEO';
                
    case 'CorrLLE'
        if length(varargin)>1, [mappedA,err ]= CorrLLE(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err]= CorrLLE(A, no_dims, varargin{1}); end
        mapping.name = 'CorrLLE';
        
    case 'CorrLNPE'
        if length(varargin)>1, [mappedA,err ]= CorrLNPE(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err]= CorrLNPE(A, no_dims, varargin{1}); end
        mapping.name = 'CorrLNPE';    
        
    case 'SafeMultiLLE'
        if length(varargin)>1, [mappedA,err ]= SafeMultiLLE(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err]= CorrLNPE(A, no_dims, varargin{1}); end
        mapping.name = 'SafeMultiLLE';    
        
    case 'MLLE'
        if length(varargin)>1, [mappedA,err]= Mlle(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err]= Mlle(A, no_dims, varargin{1}); end
        mapping.name = 'MLLE';
        
    case 'IHNE'
        if length(varargin)>1, [mappedA,err ]= IHNE(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err ] = IHNE(A, no_dims, varargin{1}); end
        mapping.name = 'IHNE';
        
    case 'RHNE'
        if length(varargin)>1, [mappedA,err ]= RHNE(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err ] = RHNE(A, no_dims, varargin{1}); end
        mapping.name = 'RHNE';
        
    case 'BHNE'
        if length(varargin)>1, [mappedA,err ]= BHNE(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err ] = BHNE(A, no_dims, varargin{1}); end
        mapping.name = 'BHNE';
        
    case 'LNPE'
        if length(varargin)>1, [mappedA,err ]= LNPE(A, no_dims, varargin{1},varargin{2});
        else [mappedA,err ] = LNPE(A, no_dims, varargin{1}); end
        mapping.name = 'LNPE';
        
    case 'LLC'
        % Compute LLC mapping
        if isempty(varargin), [mappedA,err ] = run_llc(A', no_dims, 12, 20, 200, eig_impl);
        elseif length(varargin) == 2, [mappedA,err ] = run_llc(A', no_dims, varargin{1}, 20, 200, eig_impl);
        elseif length(varargin) == 3, [mappedA,err ]= run_llc(A', no_dims, varargin{1}, varargin{2}, 200, eig_impl);
        else [mappedA,err ] = run_llc(A', no_dims, varargin{1}, varargin{2}, varargin{3}, eig_impl); end
        mappedA = mappedA';
        mapping.name = 'LLC';
        
    case 'LTSA'
        % Compute LTSA mapping
        if isempty(varargin), [mappedA,err ] = ltsa(A, no_dims, 12, eig_impl);
        else [mappedA,err ]= ltsa(A, no_dims, varargin{1}, eig_impl); end
        mapping.name = 'LTSA';
        
    otherwise
        error('Unknown dimensionality reduction technique.');
end

% JDQR makes empty figure; close it
if strcmp(eig_impl, 'JDQR')
    close(gcf);
end

% Handle PRTools dataset
if prtools == 1
    AA.data = mappedA;
    mappedA = AA;
end
