function []=solve_RE()

% load data
load teapot.mat
dataname='teapot';
teapots=X;
m=400; % 400张图像
indx=1:400;
data=teapots(:,indx(1:2:m)); % 取200个样本，每两个样本中取一个

%%%%%%%%%%%%%%%% 求重构误差 n=200        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% 柱状图可视化的实验 k=4      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% 重建误差 Table k=4~8      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% 重建误差折线图 k=2~20    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row=19;
err = zeros(3,row); % 统计各种方法的平均重建误差
for i=1:row

      K=1 + i;

     % 方法二： LLE
     [poses] = Rec_LLE(data,K);
     tempErr = sqrt(sum((poses-data).^2,1));
     err(1,i) = mean(tempErr);

     % 方法3： MLLE
     [poses] = Rec_MLLE(data,K);
     tempErr = sqrt(sum((poses-data).^2,1));
     err(2,i)=mean(tempErr);
     
     % 方法4： LNPE
     [poses] = Rec_LNPE(data,K);
     tempErr = sqrt(sum((poses-data).^2,1));
     err(3,i)=mean(tempErr);

%      % 方法4： INV_HNE
%      [poses]=Rec_IHNE(data,K);
%      tempErr = sqrt(sum((poses-data).^2,1));
%      err(3,i)=mean(tempErr);
% 
%      % 方法5： REC_HNE
%      [poses]=Rec_RHNE(data,K);
%      tempErr = sqrt(sum((poses-data).^2,1));
%      err(4,i)=mean(tempErr);
% 
%      % 方法6： BAL_HNE
%      [poses]=Rec_BHNE(data,K);
%      tempErr = sqrt(sum((poses-data).^2,1));
%      err(5,i)=mean(tempErr);

     % 打印RE
     err

end