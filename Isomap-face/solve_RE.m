function []=solve_RE()

load face_data.mat;
m=698;
data=images(:,1:m);

%%%%%%%%%%  ���ع���ƽ�����  %%%%%%%%%%%%%%%%
row=19;
err = zeros(3,row);
for i=1:row

     K = i + 1;

     % �������� LLE
     [poses] = Rec_LLE(data,K);
     tempErr = sqrt(sum((poses-data).^2,1));
     err(1,i) = mean(tempErr);

     % ����3�� MLLE
     [poses]=Rec_MLLE(data,K);
     tempErr = sqrt(sum((poses-data).^2,1));
     err(2,i)=mean(tempErr);
     
     % ����4�� LNPE
     [poses]=Rec_LNPE(data,K);
     tempErr = sqrt(sum((poses-data).^2,1));
     err(3,i)=mean(tempErr);
     
%      % ����4�� IHNE
%      [poses]=Rec_IHNE(data,K);
%      tempErr = sqrt(sum((poses-data).^2,1));
%      err(3,i)=mean(tempErr);
%      
%      % ����5�� RHNE
%      [poses]=Rec_RHNE(data,K);
%      tempErr = sqrt(sum((poses-data).^2,1));
%      err(4,i)=mean(tempErr);
%      
%      % ����6�� BHNE
%      [poses]=Rec_BHNE(data,K);
%      tempErr = sqrt(sum((poses-data).^2,1));
%      err(5,i)=mean(tempErr);

     err
     
end
