function [X, labels] = generate_data(dataname, n, d1,d2)
%GENERATE_DATA Generates an artificial dataset
%
%	[X, labels] = generate_data(dataname, n, noise)
%
% Generates an artificial dataset. Possible datasets are: 'swiss' for the Swiss roll
% dataset, 'helix' for the helix dataset, 'twinpeaks' for the twinpeaks dataset,
% '3d_clusters' for the 3D clusters dataset, and 'intersect' for the intersecting
% dataset. The variable n indicates the number of datapoints to generate 
% (default = 1000). The variable noise indicates the amount of noise that
% is added to the data (default = 0.05).
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.2b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want. However, it is appreciated if you maintain the name of the original
% author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

	if ~exist('n', 'var')
		n = 1000;
    end
    if ~exist('noise', 'var')
        noise =0;
    end

	switch dataname
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'Two Swisses'
            t = (3 * pi / 2) * (1 + 2 * rand(1,n));  
            height1 = 21 * rand(1, n); 
            height2 = 21 *(2+ rand(1,n));
           x1= [t .* cos(t); height1; t .* sin(t)];
           x2= [t .* cos(t); height2; t .* sin(t)];
           %%
            [Minvalue1,ind1]=min(x1(3,:));
            [Minvalue2,ind2]=min(x2(3,:));
            y1=ones(1,9)*x1(1,ind1)+[0.1:0.1:0.9]*(x2(1,ind2) - x1(1,ind1));
            y2=ones(1,9)*x1(2,ind1)+[0.1:0.1:0.9]*(x2(2,ind2) - x1(2,ind1));
            y3=ones(1,9)*x1(3,ind1)+[0.1:0.1:0.9]*(x2(3,ind2) - x1(3,ind1));
            %%
            X = [x1 [y1;y2;y3] x2 ]';
            %labels = [t  t]';
             labels =X(:,3);
        case 'Multi-Surfaces'
            t1 = [unifrnd(pi*5/6,pi*16/12,1,n/4)];
            t2 = [unifrnd(pi*18/12,pi*12/6,1,n/4)];
            t3=(5*pi/6)*(1+7/5*rand(1,n/2));
            a1=t1.*cos(t1); b1=t1.*sin(t1);
            c1=[unifrnd(-1,3,1,n/4)];
            a2=t2.*cos(t2); b2=t2.*sin(t2);
            c2=[unifrnd(-1,3,1,n/4)];
            a3=t3.*cos(t3); b3=t3.*sin(t3);
            c3=[unifrnd(6,10,1,n/2)];
            x1=[a1;c1;b1]; x2=[a2;c2;b2]; x3=[a3;c3;b3];
            X=[x1 x2 x3]';
            labels =X(:,3);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        case 'Two Surfaces'
            t1 = [unifrnd(pi*5/6,pi*16/12,1,n/4)];
            t2 = [unifrnd(pi*18/12,pi*12/6,1,n/4)];
            t3=(5*pi/6)*(1+7/5*rand(1,n/2));
            a1=t1.*cos(t1); b1=t1.*sin(t1);
            c1=[unifrnd(-1,3,1,n/4)];
            a2=t2.*cos(t2); b2=t2.*sin(t2);
            c2=[unifrnd(-1,3,1,n/4)];
            a3=t3.*cos(t3); b3=t3.*sin(t3);
            c3=[unifrnd(6,10,1,n/2)];
            x1=[a1;c1;b1]; x2=[a2;c2;b2]; x3=[a3;c3;b3];
            [Minvalue1,ind1]=min(x1(3,:));
            [Minvalue2,ind2]=min(x2(3,:));
            y1=ones(1,9)*x1(1,ind1)+[0.1:0.1:0.9]*(x2(1,ind2) - x1(1,ind1));
            y2=ones(1,9)*x1(2,ind1)+[0.1:0.1:0.9]*(x2(2,ind2) - x1(2,ind1));
            y3=ones(1,9)*x1(3,ind1)+[0.1:0.1:0.9]*(x2(3,ind2) - x1(3,ind1));
            X=[x1 [y1;y2;y3] x2]';
            labels =X(:,3);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'tire'
            t = pi*5*rand(1,n)/3;
            s = pi*5*rand(1,n)/3;
            X = [(3+cos(s)).*cos(t);(3+cos(s)).*sin(t);sin(s)]';
            labels = uint8(t)';
            
        case 'trefoil'
            [X] = trefoilpatch(.1,1,3,d1,d2);
            X=X';
            labels =X(:,3);
            size(labels)
       case 'S-curve'
%            t = pi*(1.5*rand(1,n/2)-1); height = 15*rand(1,n);
%            X = [[cos(t), -cos(t)]; height;[ sin(t), 2-sin(t)]]' + noise * randn(n, 3);
%            labels = uint8([cos(t), -cos(t)])';
%              angle = pi*(1.5*rand(1,n/2)-1); height = 5*rand(1,n);
%              X = [[cos(angle), -cos(angle)]; height;[ sin(angle), 2-sin(angle)]]';
%              labels =uint8([angle angle]);
           t = pi*(1.5*rand(1,n/2)-1); height = 5*rand(1,n);
            X = [[cos(t), -cos(t)]; height;[ sin(t), 2-sin(t)]]';
           labels = [t t]';
        case 'S-curve Hole'
           t = pi*(1.5*rand(1,n/2)-1); height = 5*rand(1,n);
                   for ii = 1:n/2
            if ( (t(ii) > -1.5)&(t(ii) < 0))
                if ((height(ii) > 1) & (height(ii) <3))
                    kl(ii) = 1;
                end;
            end;
        end;
        kkz = find(kl==0);
        t = t(kkz(1:n/8));
        height = height(kkz(1:n/4));
           
            X = [[cos(t), -cos(t)]; height;[ sin(t), 2-sin(t)]]';
           labels = [t t]';
           
       case 'S-curve Noise'

           t = pi*(0.5*rand(1,n/2)-1); height = 5*rand(1,n);%Ô­À´ÊÇ1.5
           X = [[cos(t), -cos(t)]+imnoise([cos(t), -cos(t)],'gaussian',0,0.1); height+imnoise(height,'gaussian',0,0.1);[ sin(t), 2-sin(t)]+imnoise([ sin(t), 2-sin(t)],'gaussian',0,0.1)]' + noise * randn(n, 3);
           labels = [t t]';

        case 'Swiss'
            t = (3 * pi / 2) * (1 + 2 * rand(1, n));  
            height = 21 * rand(1, n);
             % X = [t .* cos(t)+imnoise(t.*cos(t),'gaussian',0,0.1); height+imnoise(height,'gaussian',0,0.1); t .* sin(t)+imnoise(t.*sin(t),'gaussian',0,0.1)]' + noise * randn(n, 3);
           X = [t .* cos(t); height; t .* sin(t)]' ;
            labels = uint8(t)';
            
        case 'Swiss Hole'
        % Swiss Roll w/ hole example taken from Donoho & Grimes
        tt = (3*pi/2)*(1+2*rand(1,2*n));  
        height = 21*rand(1,2*n);
        kl = repmat(0,1,2*n);
        for ii = 1:2*n
            if ( (tt(ii) > 9)&(tt(ii) < 12))
                if ((height(ii) > 9) & (height(ii) <14))
                    kl(ii) = 1;
                end;
            end;
        end;
        kkz = find(kl==0);
        tt = tt(kkz(1:n));
        height = height(kkz(1:n));
        %X = [tt.*cos(tt)+IMNOISE(tt.*cos(tt),'gaussian',0,0.1); height+IMNOISE(height,'gaussian',0,0.1); tt.*sin(tt)+IMNOISE(tt.*sin(tt),'gaussian',0,0.1)]';     
        X = [tt.*cos(tt); height; tt.*sin(tt)]';     

        labels = uint8(tt)';
            
         case 'Swiss Connect'
        % Swiss Roll w/ hole example taken from Donoho & Grimes
        tt = (3*pi/2)*(1+2*rand(1,2*n));  
        height = 21*rand(1,2*n);
        kl = repmat(0,1,2*n);
        for ii = 1:2*n
            if ( (tt(ii) >7)&(tt(ii) < 17))
                if ((height(ii) > 9) & (height(ii) <14))
                    kl(ii) = 1;
                end;
            end;
        end;
        kkz = find(kl==0);
        tt = tt(kkz(1:n));
        height = height(kkz(1:n));
        %X = [tt.*cos(tt)+IMNOISE(tt.*cos(tt),'gaussian',0,0.1); height+IMNOISE(height,'gaussian',0,0.1); tt.*sin(tt)+IMNOISE(tt.*sin(tt),'gaussian',0,0.1)]';     
        X = [tt.*cos(tt); height; tt.*sin(tt)]';     

        labels = uint8(tt)';
        
        
        case 'Changing-Swiss'
            r = zeros(1, n);
            for i=1:n
                pass = 0;
                while ~pass
                    rr = rand(1);
                    if rand(1) > rr
                        r(i) = rr;
                        pass = 1;
                    end
                end
            end
            t = (3 * pi / 2) * (1 + 2 * r);  
            height = 21 * rand(1, n);
            X = [t .* cos(t); height; t .* sin(t)]' + noise * randn(n, 3);
            labels = uint8(t)';
            
        case 'Sphere' %by Saul & Roweis
        inc = 9/sqrt(n);   %inc = 1/4;
        [xx,yy] = meshgrid(-5:inc:5);
        rr2 = xx(:).^2 + yy(:).^2;
        [tmp ii] = sort(rr2);
        Y = [xx(ii(1:n))'; yy(ii(1:n))'];
        a = 4./(4+sum(Y.^2));
        X = [a.*Y(1,:); a.*Y(2,:); 2*(1-a)]';
        labels = X(:,3)';

            
        case 'Helix'
        	t = [1:n]'/ n;
        	t = t .^ (1.0) * 2 * pi;
			X = [(2 + cos(8 * t)) .* cos(t) (2 + cos(8 * t)) .* sin(t) sin(8 * t)]; %+ noise * randn(n, 3)
        	labels = X(:,3)';%uint8(t);
            
        case 'TwinPeaks'
            inc = 1.5 / sqrt(n);
            [xx2, yy2] = meshgrid(-1:inc:1);
            zz2 = sin(pi * xx2) .* tanh(3 * yy2);
            xy = 1 - 2 * rand(2, n);
            X = [xy; sin(pi * xy(1,:)) .* tanh(3 * xy(2,:))]' + noise * randn(n, 3);
            X(:,3) = X(:,3) * 10;
            labels = uint8(X(:,3));
            
        case '3d-Clusters'
            numClusters = 5;
            centers = 10 * rand(numClusters, 3);
            D = L2_distance(centers', centers', 1);
            minDistance = min(D(find(D > 0)));
            k = 1;
            n2 = n - (numClusters - 1) * 9;
            for i=1:numClusters
                for j=1:ceil(n2 / numClusters)
                   X(k, 1:3) = centers(i, 1:3) + (rand(1, 3) - 0.5) * minDistance / sqrt(12);
                   labels(k) = i;
                   k = k + 1;
               end
               if i < numClusters
                  for t=0.1:0.1:0.9
                       X(k, 1:3) = centers(i, 1:3) + (centers(i + 1, 1:3) - centers(i, 1:3)) * t;
                       labels(k) = 0;
                       k = k + 1;
                   end
               end
            end
            X = X + noise * randn(size(X, 1), 3);
            labels = labels';
            
        case 'Intersect'
            t = [1:n]' ./ n .* (2 * pi);
            x = cos(t);
            y = sin(t);
            X = [x x .* y rand(length(x), 1) * 5] + noise * randn(n, 3);
            labels = uint8(5 * t);
            
        case 'Intersect-2d'
            t = [1:n]' ./ n .* (2 * pi);
            x = cos(t);
            y = sin(t);
            X = [x x .* y] + noise * randn(n, 2);
            labels = uint8(5 * t);
			
		otherwise
			error('Unknown dataset name.');
	end
