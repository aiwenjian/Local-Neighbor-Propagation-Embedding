function []=displayed(D,pl,sl,faceH ,faceW,index)

fea=D; 
%numPerLine = 11; 
%ShowLine = 2; 
numPerLine = pl; 
ShowLine = sl; 
if index==1
    Y = zeros(faceW*ShowLine,faceH*numPerLine);
end
if index==2
Y = zeros(faceH*ShowLine,faceW*numPerLine);
end
for i=0:ShowLine-1 
   for j=0:numPerLine-1 
       if index==1
     Y(i*faceW+1:(i+1)*faceW,j*faceH+1:(j+1)*faceH) = reshape(fea(i*numPerLine+j+1,:),[faceW,faceH]); 
       end
       if index==2
     Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW)=reshape(fea(i*numPerLine+j+1,:),[faceH,faceW])';
       end
   end
end
if index==1
figure
imagesc(Y);
colormap(gray);
end
if index==2
   figure
   imagesc(Y);
   colormap(gray); 
end
axis off;
