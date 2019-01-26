%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **** FOR OPTION 3 ****
%
% This function processes the outputs from Jerome's Test_EWT_2D_Curvelet
% code and processes it for use in the clustering algorithm. Step by
% step, it does the following:
% - Take the k filters from Jerome's code
% - create (m times n) x k matrix of pixel values
% - calculate L2 norm in a neighborhood dependent on frequency for each
%   pixel in each filter
% 
%Input:
% - ewtc - empirical wavelet filters
% - Bw - vector of frequency values
%Output:
% - result - matrix of pixel values that can be directly given to
%   MATLAB's built-in kmeans function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [result] = ewtc2datamatrix(ewtc, Bw)

% the way the matrix is built is by initializing our kmeansmatrix to a row
% of zeros of the correct size and then using a for loop to go through each
% filter in each annulus and concatenate these pixel values at the end of
% the original matrix.

% the initial row vector:
result = zeros(size(ewtc{1}));
result = result(:)';

wedges = length(ewtc);
l = length(ewtc{1}(:));

for h = 2:(wedges)
    % we calculate the matrix for each annulus and concatenate it onto the
    % matrix we got from the center region
    k = length(ewtc{h});
    nextrows = zeros(k,l);
    entropycell = cell(1,k);
    
    Dw = zeros(length(Bw{h})+1,1);
    Dw(1) = Bw{h}(1);
    for i = 1:(length(Dw)-2)
        Dw(i+1) = (Bw{h}(i) + Bw{h}(i+1))/2;
    end
    Dw(end) = Bw{h}(end);

    for i = 1:k
    
        % turn the intensity values of the filters into local entropy values
        freq = Dw(i);
        entropycell{i} = L2nhood(ewtc{h}{i}, freq);
    
        for j = 1:l
            % put the entropy values into the k by (m x n) matrix we want
            nextrows(i,j) = entropycell{i}(j);
        end
    end
    % now we stitch it to the original matrix
    result = vertcat(result, nextrows);
end

% the built in kmeans function flips our rows/columns convention
result = result';
result = result(:, 2:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the L2 norm for each pixel in a a neighborhood
% of a size dependent on frequency
%Input:
% - I - image matrix
% - freq - frequency of the mode
%Output:
% - J - L2 norm matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = L2nhood(I, freq)
% find matrix dimensions
dims = size(I);


% Used to be 2*pi/freq, but Adina said to just go for 2/freq
T = 2*pi/freq;


% radius of the neighborhood
r = ceil(T/2);
%r = 1;
J = zeros(dims);

% some initial processing; we make everything from 0 to 255
Ipad = abs(I - mean(I(:)));
Ipad = rescale255(Ipad);



% pad the array symmetrically
% Symmetric padding messes up modes at the edge - we need to figure
% something else out

Ipad = padarray(Ipad,[r,r],'symmetric','both');

% calculate the L1 norm by subtracting out the mean value in the
% neighborhood from each pixel in the neighborhood, then sum the absolute
% values

for i = 1:dims(1)
   for j = 1:dims(2)
      K = Ipad(i:(i+2*r),j:(j+2*r));
%      K = abs(K - mean(K(:))).^2;
      J(i,j) = sum(K(:).^2)^(1/2)/((2*r+1)^2);
   end
end

end

% linearly rescales the input so that the max is 255 and the min is 0

function [rescaled, floored] = rescale255(values)
u_max = max(values(:));
u_min = min(values(:));

rescaled = 255.*(values - u_min)./(u_max - u_min);
floored = floor(rescaled);

end
