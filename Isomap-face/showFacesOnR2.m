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

%draw images on selected points
scale = 0.002;
x=zeros(64,64);
shown_images = ones(2, 100);   % 存放已经显示过的图片坐标
for p=1:size(ks,2)
    k=ks(p);
    dist = sum((poses(:, k)-shown_images).^2, 1);   % 如果当前图片与之前已经显示过的图片距离太近，则取消显示
    if min(dist) < 2e-2
        continue;
    end
    shown_images(:,p) = poses(:,k);  % 把显示过的图片放入矩阵内
    scatter(poses(1,k),poses(2,k),24,'ro');  % 把要显示的图片对应的点做标记（画红圈）
    xc=poses(1,k);
    yc=poses(2,k);
   imagesc([xc,xc+64*scale],[yc,yc+64*scale],rot90(reshape(images(:,k),64,64)'));  % 选出来的图片显示在点的左上角
   colormap gray(256);
    hold on
end
axis tight;
% axis off;
box on;

hold off
return
