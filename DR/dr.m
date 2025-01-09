function []=dr()

% This function aims to show your manifold learning method on synthetic
% manifolds, the datsets contain in the software is as follows:
%1.  'Two Swisses'
%2.  'Two Surfaces',1e-4
%3.  'Multi-Surfaces'
%4.  'trefoil'1e-3 %d2=1(d1=0.5 ,532个点，d1=0.4,665),d2=3(d1=0.1,1131个点；d1=0.2 ,567个点；d1=0.3,378个点,d1=0.8,144个点)(最后两个参数可调)
%5.  'S-curve',1e-4;
%6.  'Swiss',1e-6
%7.  'Swiss Hole',1e-4
%8.  'Changing-Swiss',1e-4
%9. 'Sphere',1e-2
%10. 'Helix',1e-3
%11. 'Swiss Connect'
%12. '3d-Clusters'
%13. 'Intersect'

dataset = {'S-curve','Swiss','Swiss Hole','Changing-Swiss','Sphere','Helix','Swiss Connect','TwinPeaks','3D-Clusters','Intersect'};
tol = [1e-3,1e-4,1e-4,1e-3,1e-3,1e-4,1e-4,1e-4,1e-4,1e-2,1e-2,1e-4,1e-3,1e-4,1e-3];

numpoints = 600;            % 生成数据点的个数

% 选择需要对比的方法
%alg = {'LLE', 'HLLE', 'LLC', 'MLLE', 'LTSA', 'LNPE', 'CorrLLE'};  
alg = {'LLE', 'SafeMultiLLE', 'CorrLNPE', 'CorrLLE'}; 
col = length(alg);    %方法个数 6
row = 5;                    % 画图的行数
no_dims =2;                 % 低维维度

for p = 1 : 6  % 选择数据集
    dataname=dataset{p};
    if strcmp(dataname,'trefoil')
        [X, labels] = generate_data(dataname,numpoints,0.8,3);  %生成数据集 数据赋给X

    else
        [X, labels] = generate_data(dataname,numpoints);

    end 
%     load('2.mat')
%     X=Y(:,1:3);
%     labels = Y(:,4);
    
    % 画图。。。。。。。。。。。。 三维 原式数据集
    figure,
    scatter3(X(:,1), X(:,2), X(:,3), 23, labels,'.');  
    axis tight;
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
    grid off;
    box on;
    zlabel(dataname,'FontSize',15);
    figure
    % 画图。。。。。。。。。。。。三维  原式数据集
    
    for i=1:row
        k=3+i;              % 设置k的起始大小
        for j=1:col
            subplot(row,col,(i-1)*col+j)
            [mappedX]= compute_mapping(X, alg{j}, no_dims, k,tol(p));  %调用降维方法
            scatter(mappedX(:,1), mappedX(:,2), 23, labels,'.');
            axis tight;
            set(gca,'xtick',[],'ytick',[]);
            box on;
            if j==1
                ylabel(strcat('k=',num2str(k)),'FontSize',15);
            end
            if i==row
                xlabel(alg{j},'FontSize',15);
            end
        end
    end
end