%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a script that performs cartoon+texture decomposition and then
%cartoon image segmentation using multiphase and local multiphase 
%segmentation. This runs on all of the test images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the images
load('images.mat');

%format the figure
iptsetpref('ImshowBorder','tight');
for file = 1:16
    %obtain an image
    im_name = strcat('image', num2str(file));
    
    %perform cartoon+texture decomposition
    [u, ~] = CartoonTexture_nonlinear(eval(im_name), 3);
    
    %perform multilevel image threshold by Otsu's method
    thresh = multithresh(u, 3);
    thresh_result = imquantize(u, thresh);

    %perform segmentation using k-means
    k_mean_result = imsegkmeans(uint8(255*mat2gray(u)),4);
    
    %perform multiphase segmentation
    fprintf('Performing multiphase segmentation on image %d.\n', file)
    pam1 = cartoon_parameters{file,1};
    m_result = multiphase(u, pam1(1), pam1(2), pam1(3));
    
    %perform local multiphase segmentation
    fprintf('Performing local multiphase segmentation on image %d.\n', file)
    pam2 = cartoon_parameters{file,2};
    lm_result = localmultiphase(u, pam2(1), pam2(2), pam2(3), pam2(4));
    
    %output the images
    figure; 
	subplot(2,3,1); imagesc(u); title('cartoon'); axis off; axis square; colormap gray;
	subplot(2,3,2); imagesc(thresh_result-1); title('Otsu Method'); axis off; axis square; colormap gray;
	subplot(2,3,3); imagesc(k_mean_result); title('k means'); axis off; axis square; colormap gray;
	subplot(2,3,4); imagesc(m_result); title('multiphase'); axis off; axis square; colormap gray;
	subplot(2,3,5); imagesc(lm_result); title('local multiphase'); axis off; axis square; colormap gray;
    
end