function poltdr(dataname,numpoints,alg1,alg2,alg3,alg4,k)
    if strcmp(dataname,'trefoil')
    [X, labels] = generate_data(dataname,numpoints,0.2,3);
    else
       [X, labels] = generate_data(dataname,numpoints); 
    end

    subplot(1,5,1)
    scatter3(X(:,1), X(:,2), X(:,3), 23, labels,'.'); drawnow%title('Original dataset'), 

    title(strcat('Origional-',dataname,',','k=',num2str(k)));
	no_dims =2; %round(intrinsic_dim(X, 'MLE'));
	disp(['MLE estimate of intrinsic dimensionality:' num2str(no_dims)]);
	mappedX = compute_mapping(X, alg1, no_dims, k);
   

   subplot(1,5,2)
   scatter(mappedX(:,1), mappedX(:,2), 23, labels,'.');
   xlabel('LLC');
 
 
    %compare option
    if strcmp(alg2, 'NN')
        disp('No Compare option')
    else
    mappedX1 = compute_mapping(X, alg2, no_dims,k);   
    
    subplot(1,5,3)
       scatter(mappedX1(:,1), mappedX1(:,2),23,labels,'.'); %title('Result of dimensionality reduction opt'), drawnow
       xlabel('LTSA');
     end

     mappedX2 = compute_mapping(X, alg3, no_dims,k);

       subplot(1,5,4)    
       scatter(mappedX2(:,1), mappedX2(:,2), 23, labels,'.'); %title('Result of dimensionality reduction opt'), drawnow
       xlabel('LLE');

     mappedX3 = compute_mapping(X, alg4, no_dims,k);
     
       subplot(1,5,5)
       scatter(mappedX3(:,1), mappedX3(:,2), 23, labels,'.'); %title('Result of dimensionality reduction opt'), drawnow
       xlabel('STLLE');
       return