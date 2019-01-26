%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a script that performs cartoon+texture decomposition and then
%texture image segmentation by clustering the L2 energy matrix of each
%image of the filter bank computed from Curvelet Empirical Wavelet 
%Transform. This runs on a specified image from the set of test images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load images
load('images.mat');

%specify image number here
image = 1;

%set up parameters for EWTC
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

%print message
fprintf('Performing texture segmentation on image %d.\n', image);

%obtain an image
im_name = strcat('image',num2str(image));

%obtain parameters
pam = texture_parameters{image};

%perform cartoon+texture decomposition
[~,v] = CartoonTexture_nonlinear(eval(im_name),3);

%Perform EWT
[ewtc, mfb, Bw,Bt] = modified_EWT2D_Curvelet(v, params, pam(1));

%obtain data matrix
result = ewtc2datamatrix(ewtc, Bw);

%perform clustering
[s,k] = texture_cluster(result, pam(2), pam(3));

%perform texture segmentation by gabor filter
[X, g] = gabor_cluster(v, pam(2));

%output the images
figure; 
subplot(2,2,1); imagesc(v); axis off; axis square; colormap gray; title('texture');
subplot(2,2,2); imagesc(reshape(g, 256, 256)); axis off; axis square; colormap gray; title("Gabor Filter");
subplot(2,2,3); imagesc(reshape(k,256,256)); axis off; axis square; colormap gray; title('k-means');
subplot(2,2,4); imagesc(reshape(s,256,256)); axis off; axis square; colormap gray; title('MBO');