%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs Empirical Wavelet Transform with Curvelet (EWTC) 
%onto the image. After the filter bank from EWTC is obtained, a data
%matrix is created based on it by calculating the entropy values of each
%pixel on each image of the filter bank.
%
%Input:
% image - image that will be converted to a data matrix
%Output:
% result - data matrix of the image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = obtain_datamatrix(image)

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


%perform EWTC onto the image
[ewtc, ~, Bw] = modified_EWT2D_Curvelet(image, params);

%obtain the data matrix
result = ewtc2datamatrix(ewtc, Bw);

end