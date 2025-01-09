function [Y, T] = BHNE(data, K, d)

    X = data;  %ÿ����һ�����ݵ�
    [D, N] = size(X);
    X2 = sum(X.^2,1);
    tol=1e-5;
    distance = repmat(X2, N, 1) + repmat(X2', 1, N) - 2 * X' * X;
    [sorted, index] = sort(distance);
    neighborhood = index(2:(1+K), :);
    
    %���ó�����
    EPOCH = 2;  % �����ִ�
    
    %++++++++++++++++++++++++++++  ���LLEȨֵ  +++++++++++++++++++++++++++
    fprintf(1,'-->Solving for reconstruction weights.\n');
    W = zeros(K,N);
    for ii=1:N
       z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K);     % shift ith pt to origin
       C = z'*z;                                            % local covariance
       C = C + eye(K,K)*tol*trace(C);                       % regularlization (K>D)     
       W(:,ii) = C\ones(K,1);                               % solve Cw=1
       W(:,ii) = W(:,ii)/sum(W(:,ii));                      % enforce sum(w)=1
    end

    % ++++++++++++++++++++++++++ dirAMLȨֵ������  ++++++++++++++++++++++
    fprintf(1,'-->Solving for 2nd reconstruction weights.\n');
    T=zeros(N, N);
    tol=1e-6;
    for ii=1:N
        R=zeros(N, K);
%         z = X(:, neighborhood(:, neighborhood(:, ii))) - repmat(X(:,ii), 1, K*K);    % shift ith pt to origin
%         C = z' * z;                                                                  % local covariance,shape(K^2, k^2)
%         C = C + eye(K*K, K*K) * tol * trace(C);                                      % regularlization (K>D)
%         WW = C \ ones(K*K, 1);                                                       % solve Cww=1
%         WW = WW / sum(WW);                                                           % enforce sum(ww)=1
        
    %+++++++++++++++++++++++++++++++ ��ÿ������ؽ����е��� ++++++++++++++++++++++++++++++
        
        %�ø��ºõ�ȥ������һ��
%         for epoch = 1:EPOCH  
%             for j = 1:K                                                              %j��Xi�ĵ�j�����ڵ�
%                 recX_i = X(:, neighborhood(:, neighborhood(:, ii))) * WW;            % recX_i������ֱ��Ȩֵ�����������ؽ�X_i
%                 ori_WW = WW((j-1)*K+1:j*K);                                          %������Ҫ���µ�WWȨֵ
%                 objX_ij = X(:, neighborhood(:, neighborhood(j, ii))) * ori_WW;
%                 X_ij = X(:,ii) - (recX_i - objX_ij);
%                 X_ij_ =  W(j, ii) * X(:, neighborhood(:, neighborhood(j, ii))) - repmat(X_ij, 1, K);
%                 CC = X_ij_' * X_ij_;
%                 CC = CC + eye(K, K) * tol * trace(CC);
%                 W_ij = CC \ ones(K, 1);
%                 W_ij = W_ij / sum(W_ij);                                             %������ĵڶ���Ȩֵ
%                 WW((j-1)*K+1:j*K) = W_ij * W(j, ii);
%             end
%         end

         %��LLEΪԭ����Ϊ��һ�ֵĵ����������ӵڶ��ֿ�ʼ����ԭʼ����
         WW = zeros(K*K, 1);
         tmp = [];
         for epoch = 1:EPOCH
            for j = 1:K                                                              %j��Xi�ĵ�j�����ڵ�
                if epoch == 1                                                        %��һ�ֻ���LLE
                    W_i = W(:,ii) ; W_i(j)=0;
                    recX_i = X(:, neighborhood(:, ii)) * W_i;
                    X_ij = X(:,ii) - recX_i;                                           %����Ҫ���µĵ�ֵ
                else
                    recX_i = X(:, neighborhood(:, neighborhood(:, ii))) * WW;            % recX_i������ֱ��Ȩֵ�����������ؽ�X_i
                    ori_WW = WW((j-1)*K+1:j*K);                                          %������Ҫ���µ�WWȨֵ
                    objX_ij = X(:, neighborhood(:, neighborhood(j, ii))) * ori_WW;
                    X_ij = X(:,ii) - (recX_i - objX_ij);
                end  
                X_ij_ = W(j, ii) * X(:, neighborhood(:, neighborhood(j, ii))) - repmat(X_ij, 1, K);
                CC = X_ij_' * X_ij_;
                CC = CC + eye(K, K) * tol * trace(CC);
                W_ij = CC \ ones(K, 1);
                W_ij = W_ij / sum(W_ij);                                                %������ĵڶ���Ȩֵ
                WW((j-1)*K+1:j*K) = W_ij * W(j, ii);
            end
            recX_i = X(:, neighborhood(:, neighborhood(:, ii))) * WW;
            tmp(epoch) = norm(X(:, ii) - recX_i);
            
            
         end
         tmp;

%         %�ó�ʼ��ȥ������һ��
%         newWW = zeros(K*K, 1);
%         for epoch = 1:EPOCH
%             for j = 1:K                                                              %j��Xi�ĵ�j�����ڵ�
%                 recX_i = X(:, neighborhood(:, neighborhood(:, ii))) * WW;            % recX_i������ֱ��Ȩֵ�����������ؽ�X_i
%                 ori_WW = WW((j-1)*K+1:j*K);                                          %������Ҫ���µ�WWȨֵ
%                 objX_ij = X(:, neighborhood(:, neighborhood(j, ii))) * ori_WW;
%                 X_ij_ = recX_i - objX_ij - repmat(X(:,ii), 1, K);
%                 CC = X_ij_' * X_ij_;
%                 CC = CC + eye(K, K) * tol * trace(CC);
%                 W_ij = CC \ ones(K, 1);
%                 W_ij = W_ij / sum(W_ij);                                             %������ĵڶ���Ȩֵ
%                 newWW((j-1)*K+1:j*K) = W_ij * W(j, ii);
%             end
%             WW = newWW;
%         end
        
        reshapeWW = reshape(WW, K, K) ./ repmat(W(:, ii)', K, 1);
        for jj = 1:K
            R(neighborhood(:,neighborhood(jj,ii)),jj)=reshapeWW(:, jj);
        end
        T(:,ii)=R*W(:,ii);
    end

M2 = sparse(1:N,1:N,ones(1,N),N,N,4*K*N); 
for ii=1:N
   w = W(:,ii);
   jj = neighborhood(:,ii);
   M2(ii,jj) = M2(ii,jj) - w';
   M2(jj,ii) = M2(jj,ii) - w;
   M2(jj,jj) = M2(jj,jj) + w*w';
end;

%%+++++++++++++++  Low Dimensional Embedding +++++++++++++++++++++++++++
M=(eye(N,N)-T)*(eye(N,N)-T)'+full(M2);
%%%%%%%%%%%%%%% CALCULATION OF EMBEDDING
% options.disp = 0; options.isreal = 1; options.issym = 1; 
% [Y,eigenvals] = eigs(M,d+1,0,options);
if condest(M)==inf  %strcmp(eig_impl, 'JDQR')
    disp('SINGULAR MATRIX');
    options.Disp = 0;
    options.LSolver ='bicgstab' ;%'gmres'
    [Y, eigenvals] = jdqr(M, d + 1, tol, options);      % only need bottom (no_dims + 1) eigenvectors
else

options.disp = 0; options.isreal = 1; options.issym = 1; 
[Y,eigenvals] = eigs(M,d+1,0,options);

end

[eigenvals, ind] = sort(diag(eigenvals), 'ascend');
Y = Y(:,ind(2:d+1))*sqrt(N);                                                          % bottom evect is [1,1,1,1...] with eval 0
Y = Y';
fprintf(1,'Done.\n');