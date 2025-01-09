function drawfig()

load teapot.mat
teapots=255*X;
d=2;
m=400;
indx=1:400;
teapots=teapots(:,indx(1:2:m));
row=1;
col=5;
ks=1:4:size(teapots,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  降维实验  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:row
    
    K=4;

    % 1. LLE
    [poses]=lle(teapots,K,d); 
    subplot(row,col,(i-1)*col+1)
    showFacesOnR2(teapots,poses,ks);

    % % %  % 2. HLLE 未使用，跑不动
    % % %   [poses] = hlle(teapots', d, K, 'eig_impl');
    % % %  subplot(row,col,(i-1)*col+2)
    % % %  showFacesOnR2(teapots,poses,ks);

    %  2. LLC
    [poses] = run_llc(teapots, d, K, 1, 200,'eig_impl');
    subplot(row,col,(i-1)*col+2)
    showFacesOnR2(teapots,poses,ks);

    % 3. LTSA
    [poses,err] = ltsa(teapots, d, K);
    subplot(row,col,(i-1)*col+3)
    showFacesOnR2(teapots,poses,ks);

    %  4. MLLE
    [poses,err] = Mlle(teapots, d, K);
    subplot(row,col,(i-1)*col+4)
    showFacesOnR2(teapots,poses,ks);

    % 5. LNPE
    [poses,err] = LNPE(teapots, d, K);
    subplot(row,col,(i-1)*col+5)
    showFacesOnR2(teapots,poses,ks);

end