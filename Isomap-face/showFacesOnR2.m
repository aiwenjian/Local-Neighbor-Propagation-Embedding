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
shown_images = ones(2, 100);   % ����Ѿ���ʾ����ͼƬ����
for p=1:size(ks,2)
    k=ks(p);
    dist = sum((poses(:, k)-shown_images).^2, 1);   % �����ǰͼƬ��֮ǰ�Ѿ���ʾ����ͼƬ����̫������ȡ����ʾ
    if min(dist) < 2e-2
        continue;
    end
    shown_images(:,p) = poses(:,k);  % ����ʾ����ͼƬ���������
    scatter(poses(1,k),poses(2,k),24,'ro');  % ��Ҫ��ʾ��ͼƬ��Ӧ�ĵ�����ǣ�����Ȧ��
    xc=poses(1,k);
    yc=poses(2,k);
   imagesc([xc,xc+64*scale],[yc,yc+64*scale],rot90(reshape(images(:,k),64,64)'));  % ѡ������ͼƬ��ʾ�ڵ�����Ͻ�
   colormap gray(256);
    hold on
end
axis tight;
% axis off;
box on;

hold off
return
