%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the code for the segmentation framework. Running this code will
%show the original image, the cartoon image and its segmentation, and the
%texture image and its segmentation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%Load the Image Here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify the image here
image = eval(strcat('image', num2str(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%Cartoon Segmentation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%perform cartoon+texture decomposition
fprintf('Performing Cartoon+Texture Decomposition.\n');
[u, v] = CartoonTexture_nonlinear(image, 3);
fprintf('Cartoon+Texture Decomposition Completed.\n\n');

%specify parameters for local multiphase segmentation
lambda = 10; %weighing parameter for the fidelity term
beta = 10; %weighing parameter for the intensity difference term
nu = 255^2*10^(-3); %weighing parameter for the regularization term
dt = 3.20; %time step

%perform local multiphase segmentation
fprintf('Performing Cartoon Segmentation.\n');
cartoon_result = localmultiphase(u, lambda, beta, nu, dt);
fprintf('Cartoon Segmentation Completed.\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%Texture Segmentation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set up parameters for empirical wavelet transform (curvelet)
params.option=3; 
params.typeDetect='otsu';
params.log=0;
params.detect = 'scalespace';
params.globtrend = 'none';
params.reg = 'none';
params.completion = 0;
params.curvpreproc = 'none';
params.curvreg = 'none';
params.curvdegree = 0;
params.curvlengthFilter = 0;
params.curvsigmaFilter = 0;
params.curvmethod = 'scalespace';

%specify parameters for texture clustering
tau = 0.92; %thresholding parameter (percentile)
k = 5; %number of clusters
dt = 0.03; %time step

%Run empirical wavelet transform (curvelet)
fprintf('Performing Texture Segmentation.\n');
[ewtc, mfb, Bw,Bt] = modified_EWT2D_Curvelet(v, params, tau);

%obtain data matrix
result = ewtc2datamatrix(ewtc, Bw);

%perform clustering
s = texture_cluster(result, k, dt);
texture_result = reshape(s, size(image,1), size(image,2));
fprintf('Texture Segmentation Completed.\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output the segmentation results
figure;
subplot(2,3,1); imagesc(image); axis off; axis square; colormap gray; title('Original');
subplot(2,3,2); imagesc(u); axis off; axis square; colormap gray; title('Cartoon');
subplot(2,3,3); imagesc(cartoon_result); axis off; axis square; colormap gray; title('Cartoon Segmentation');
subplot(2,3,5); imagesc(v); axis off; axis square; colormap gray; title('Texture');
subplot(2,3,6); imagesc(texture_result); axis off; axis square; colormap gray; title('Texture Segmentation');
