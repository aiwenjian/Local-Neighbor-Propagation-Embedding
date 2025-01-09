function drawfig()

load face_data.mat
d=2;
m=698;  %数据点数
data=images(:,1:m);
% k = [ 5, 6, 8];
row=4;
col=5;  
ks=1:4:m;

for i=1:row
    
    K=4 + i;

    % 1. LLE
    [poses]=lle(data,K,d); 
    subplot(row,col,(i-1)*col+1)
    showFacesOnR2(images,poses,ks);

    % 2. LLC
    [poses] = run_llc(data, d, K, 1, 200,'eig_impl');
    subplot(row,col,(i-1)*col+2)
    showFacesOnR2(images,poses,ks);

    % 3. LTSA
    [poses,err] = ltsa(data, d, K);
    subplot(row,col,(i-1)*col+3)
    showFacesOnR2(images,poses,ks);

    % 4. MLLE
    [poses,err] = Mlle(data, d, K);
    subplot(row,col,(i-1)*col+4)
    showFacesOnR2(images,poses,ks);
    
    % 5. LNPE
    [poses,err] = LNPE(data, d, K);
    subplot(row,col,(i-1)*col+5)
    showFacesOnR2(images,poses,ks);

%     % 5. IHNE
%     [poses] = IHNE(data,K,d);
%     subplot(row,col,(i-1)*col+5)
%     showFacesOnR2(images,poses,ks);
% 
%     % 6. RHNE
%     [poses, ~] = RHNE(data,K,d);
%     subplot(row,col,(i-1)*col+6)
%     showFacesOnR2(images,poses,ks);
% 
%     % 7. BHNE
%     [poses, ~] = BHNE(data,K,d);
%     subplot(row,col,(i-1)*col+7)
%     showFacesOnR2(images,poses,ks);
 
end