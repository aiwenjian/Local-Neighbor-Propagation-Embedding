function []=hist_experiment()

load face_data.mat;
m=698;
data=images(:,1:m);
row=1;

%%%%%%%%%%%%%%%%%%%%%%%  hist��Ӧ��ͼ������ͼ���ؽ����  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  LNPE�ؽ�Ч��һ�㣬���Ǽ���Ŀ죬��LLEһ������
for i=1:row

    K = 5;
    M = zeros(4, 4096);
    m_error = zeros(1, 3);
    idx = 568;              % ѡ��ǰ��������

    M(1, :) = uint8([255*data(:,idx)']);

    [poses] = Rec_LLE(data,K);
    m_error(1) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
    M(2, :) = uint8([255*poses(:,idx)']);

    [poses]=Rec_MLLE(data,K);
    m_error(2) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
    M(3, :) = uint8([255*poses(:,idx)']);
    
    [poses]=Rec_LNPE(data,K);
    m_error(3) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
    M(4, :) = uint8([255*poses(:,idx)']);

%     [poses]=Rec_IHNE(data,K);
%     m_error(3) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
%     M(4, :) = uint8([255*poses(:,idx)']);
% 
%     [poses]=Rec_RHNE(data,K);
%     m_error(4) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
%     M(5, :) = uint8([255*poses(:,idx)']);
% 
%     [poses]=Rec_BHNE(data,K);
%     m_error(5) = sqrt(sum((poses(:,idx)-data(:,idx)).^2));
%     M(6, :) = uint8([255*poses(:,idx)']);

    displayed(M, 4, 1, 64, 64, 1);
    m_error

end
