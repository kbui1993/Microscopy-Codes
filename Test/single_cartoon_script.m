%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a script that performs cartoon+texture decomposition and then
%cartoon image segmentation using multiphase and local multiphase 
%segmentation. This runs on a specified image of the set of test images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the images
load('images.mat');

%specify the image here
image = 2;

%obtain the image
im_name = strcat('image', num2str(image));

%perform cartoon+texture decomposition
[u, ~] = CartoonTexture_nonlinear(eval(im_name), 3);

%perform multiphase segmentation
fprintf('Performing multiphase segmentation on image %d.\n', image)
pam1 = cartoon_parameters{image,1};
m_result = multiphase(u, pam1(1), pam1(2), pam1(3));

%perform local multiphase segmentation
fprintf('Performing local multiphase segmentation on image %d.\n', image)
pam2 = cartoon_parameters{image,2};
lm_result = localmultiphase(u, pam2(1), pam2(2), pam2(3), pam2(4));

%format the figure
iptsetpref('ImshowBorder','tight');

%output the images
figure; 
subplot(1,3,1); imagesc(u); title('cartoon'); axis off; axis square; colormap gray;
subplot(1,3,2); imagesc(m_result); title('multiphase'); axis off; axis square; colormap gray;
subplot(1,3,3); imagesc(lm_result); title('local multiphase'); axis off; axis square; colormap gray;




