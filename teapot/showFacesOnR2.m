function showFacesOnR2(images,poses,ks);
%normalize into 1:1
poses(1,:)=poses(1,:)/range(poses(1,:));
poses(2,:)=poses(2,:)/range(poses(2,:));
%draw all points
scatter(poses(1,:),poses(2,:),12,'ko','filled');
% xlabel('left-right pose');
% ylabel('up-down pose'); 
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
hold on
%draw selected points
scatter(poses(1,ks),poses(2,ks),24,'ro');
hold on
%draw images on selected points
scale = 0.001;

for p=1:size(ks,2)
    k=ks(p);
    xc=poses(1,k);
    yc=poses(2,k);
   %imagesc([xc,xc+64*scale],[yc,yc+64*scale],reshape(images(:,k),76,101,3)/255);
   imagesc([xc,xc+64*scale],[yc,yc-64*scale],reshape(images(:,k),76,101,3)/255);
   colormap gray(256);
    hold on
end
axis tight;
% axis off
box on;
hold off
return
