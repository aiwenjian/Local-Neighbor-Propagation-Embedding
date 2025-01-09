function []=displayed(D,pl,sl,faceH ,faceW,index)
fea=D;
 
%numPerLine = 11; 
%ShowLine = 2; 
numPerLine = pl; 
ShowLine = sl; 
if index==1
    Y = zeros(faceW*ShowLine,faceH*numPerLine,3); 
end
if index==2
Y = zeros(faceH*ShowLine,faceW*numPerLine,3);
end
for i=0:ShowLine-1 
   for j=0:numPerLine-1 
       if index==1
     Temp1=reshape(fea(i*numPerLine+j+1,:),[faceW,faceH,3]);
     Y(i*faceW+1:(i+1)*faceW,j*faceH+1:(j+1)*faceH,1) = Temp1(:,:,1); 
     Y(i*faceW+1:(i+1)*faceW,j*faceH+1:(j+1)*faceH,2) = Temp1(:,:,2); 
     Y(i*faceW+1:(i+1)*faceW,j*faceH+1:(j+1)*faceH,3) = Temp1(:,:,3);  
       end
       if index==2
     Temp2=reshape(fea(i*numPerLine+j+1,:),[faceH,faceW,3]);
     Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW,1)=[Temp2(:,:,1)]'; 
     Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW,2)=[Temp2(:,:,2)]'; 
     Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW,3)=[Temp2(:,:,3)]'; 
       end
   end
end
if index==1
figure
imagesc(Y/255);
set(gca,'xtick',[],'ytick',[]);
end
if index==2
   figure
   imagesc(Y/255);
   set(gca,'xtick',[],'ytick',[]);
end
axis off
