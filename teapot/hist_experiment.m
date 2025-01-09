function []=hist_experiment()

% load data
load teapot.mat
dataname='teapot';
teapots=X;
m=400; % 400��ͼ��
indx=1:400;
data=teapots(:,indx(1:2:m));

%%%%%%%%%%%%%%% ���ֱ��ͼ���ӻ���ʵ�� k = 4, n = 200 %%%%%%%%%%%%%%%%
row=1;
for i=1:row
    
    K=4;

    M = zeros(4, 23028); % ����ؽ���ĸ�ά����
    m_error = zeros(1, 3); % ͳ�Ƶ�ǰidxͼ����ؽ����
    idx = 141; % ѡ��ǰ���ĸ�index��ͼ��

    % ԭͼ
    M(1, :) = uint8([255*data(:,idx)']); 

    % LLE
    [poses] = Rec_LLE(data,K);
%     %%%%%%%%% LLE���ؽ����������򣬴�ӡerror�Ͷ�Ӧ��index   begin %%%%%%%%%%
%     [LLE_error, index] = sort(sqrt(sum((poses(:,:)-data(:,:)).^2, 1)));
%     LLE_error
%     index
%     %%%%%%%%% LLE���ؽ����������򣬴�ӡerror�Ͷ�Ӧ��index   end   %%%%%%%%%%
    m_error(1) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
    M(2, :) = uint8([255*poses(:,idx)']);

    % MLLE
    [poses]=Rec_MLLE(data,K);
    m_error(2) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
    M(3, :) = uint8([255*poses(:,idx)']);
    
    % LNPE
    [poses]=Rec_LNPE(data,K);
    m_error(3) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
    M(4, :) = uint8([255*poses(:,idx)']);

%     % IHNE
%     [poses]=Rec_IHNE(data,K);
%     m_error(3) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
%     M(4, :) = uint8([255*poses(:,idx)']);
% 
%     % RHNE
%     [poses]=Rec_RHNE(data,K);
%     m_error(4) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
%     M(5, :) = uint8([255*poses(:,idx)']);
% 
%     % BHNE
%     [poses]=Rec_BHNE(data,K);
%     m_error(5) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
%     M(6, :) = uint8([255*poses(:,idx)']);

    displayed(M, 4, 1, 101,76, 1);
    m_error % ��ӡerror
     
end
