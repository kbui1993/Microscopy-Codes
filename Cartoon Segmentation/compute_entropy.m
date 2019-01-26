%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function computes the entropy measure of the segmentation. The
%entropy measure is provided in 
%http://cs.slu.edu/~fritts/papers/spie04_entropy_segeval.pdf.
%
%Input:
%   -image: image being segmented
%   -cluster: segmentation
%
%output:
%   -s_H_r: expected region entropy
%   -S_H_l: layout entropy
%   -entropy: total entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s_H_r, s_H_l, entropy] = compute_entropy(image, cluster)
%reshape image and segmentation as vectors
npixel = size(image,1)*size(image,2);
image = reshape(uint8(255*mat2gray(image)), npixel,1);
cluster = reshape(cluster+1, npixel,1);

%initialize variables to assist in computation of entropy
H_r = zeros(4,1);
H_l = zeros(4,1);

%compute entropy for each region
for i = 1:4
    %obtain the region
    ind = (cluster == i);
    region = image(ind);
    
    %compute expected region entropy for a region
    L = hist(double(region), unique(double(region)));
    H_R = -sum((L/sum(ind)).*log(L/sum(ind)));
    H_r(i) = (sum(ind)/npixel)*H_R;
    
    %compute layout entropy for a region
    H_l(i) = -(sum(ind)/npixel)*log(sum(ind)/npixel);
end

%compute total expected region entropy and layout entropy
s_H_r = sum(H_r);
s_H_l = sum(H_l);

%add up the two entropy
entropy = s_H_r + s_H_l;
end