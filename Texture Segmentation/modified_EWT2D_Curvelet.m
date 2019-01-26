%==========================================================================
% function [ewtc,mfb,Bw,Bt] = modified_EWT2D_Curvelet(f,params, tau)
%
% This function performs the Empirical Curvelet Transform. The Fourier 
% boundaries and angles are detected using the Pseudo-Polar FFT. 
%
% TO RUN THIS FUNCTION YOU NEED TO HAVE THE MATLAB POLARLAB TOOLBOX OF 
% MICHAEL ELAD: http://www.cs.technion.ac.il/~elad/Various/PolarLab.zip
%
% Input:
%   -f: input image
%   -params: structure containing the following parameters:
%       -params.log: 0 or 1 to indicate if we want to work with
%                    the log of the ff
%       -params.typeDetect: (for scalespace method only) 'otsu',
%                           'halfnormal','empiricallaw','mean','kmeans'
%       -params.option: 1 = scales and angles are independants
%                       2 = scales first then angles per scales
%                       3 = angles first then scales per angle
%   -tau: thresholding parameter between [0,1]. Set the lowest tau
%          percentage of entries with respect to the absolute values of 
%          Fourier coefficients to zero.
%
% Output:
%   -ewtc: cell containing each filtered output subband (ewtc{1} is the 
%   lowpass subband the next ewtc{s}{t} are the bandpass filtered images, s
%   corresponds to the scales and t to the direction)
%   -mfb: cell containing the set of empirical filters in the Fourier 
%   domain (the indexation is the same as ewtc above)
%   -Bw: list of the detected scale boundaries
%   -Bt: list of the detected angle boundaries
%   -BwN: list of the NON normalized detected scale boundaries
%   -BtN: list of the NON normalized detected angle boundaries
%
% This is a modification to the code by
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2017
%==========================================================================

function [ewtc,mfb,Bw,Bt,BtN,BwN] = modified_EWT2D_Curvelet(f,params, tau)

W=size(f,2);
H=size(f,1);


% Pseudo Polar FFT of f
PseudoFFT=PPFFT(f);

% Apply threshold to PseudoFFT
percentile = prctile(abs(PseudoFFT(:)), tau*100); 
PseudoFFT(abs(PseudoFFT) <=percentile) = 0;

%% Modified Version of Option 3 - angles first and then scales per angle
        
% Compute the mean spectrum with respect to the angle to find the first scale
meanppfft = fftshift( sum( abs(PseudoFFT) ,2));
mid_index = round( length(meanppfft)/2 );  


% Detect the first scale boundary
boundariesW = EWT_Boundaries_Detect(meanppfft(mid_index:end),params);
BwN1 = boundariesW(1);
Bw1 = boundariesW(1)*pi/mid_index;

% Compute the mean spectrum with respect to the magnitude frequency to find
% the angles
start_of_textures = boundariesW(1) + floor( size(PseudoFFT,1)/2 );

meanppfft=sum( abs( PseudoFFT(start_of_textures:end,:)),1);

%Detect the boundaries
boundaries = EWT_Angles_Detect(meanppfft', params);

BtN = (boundaries-1);
Bt = (boundaries-1)*pi/length(meanppfft)-3*pi/4;

Bw=cell(length(boundaries)+1,1);
Bw{1}=Bw1;
BwN=cell(length(boundaries)+1,1);
BwN{1}=BwN1;

% we detect the scales per angular sector
for t=1:length(boundaries)-1
    %average spectrum on the given angular range
    meanppfft=sum(abs(PseudoFFT(start_of_textures:end,boundaries(t):boundaries(t+1))),2);
    % Detect the boundaries

    bounds = EWT_Boundaries_Detect(meanppfft,params);
    BwN{t+1} = boundariesW(1)+bounds;
    Bw{t+1} = (boundariesW(1)+bounds)*pi/mid_index;
end
%last one
%average spectrum on the given angular range
meanppfft=...
    sum(abs(PseudoFFT(start_of_textures:end,boundaries(end):end)),2)+...
    sum(abs(PseudoFFT(start_of_textures:end,1:boundaries(1))),2);

% Detect the boundaries

bounds = EWT_Boundaries_Detect(meanppfft,params);
BwN{end} = boundariesW(1)+bounds;
Bw{end} = (boundariesW(1)+bounds)*pi/mid_index;


for i=2:length(Bw)
    if isempty(Bw{i})
        Bw{i}(1)=[Bw{1}];
    end

end
%% NEW CODE 
% Polar matrix where modes are located
PolarMatrix = fftshift(abs(PseudoFFT));

% temporary scale lengths
Scale_lengths = Bw;


% Go to every sector and find the L1 (sum everything and divide by
% the area) First we must find the max L1 to then set a threshold

% initialize max
max_L1= 0;


% Compute L1 and find max
for a=1:(length(BtN)-1)
    for s=2:(length(BwN)-1)
        % Go from the outermost scale towards the center
        for t=length(BwN{s}): -1: 2

            width = BwN{s}(t) -BwN{s}(t-1);
            height = BtN(a+1)-BtN(a);

            area_sector= width*height;

            sector_sum = sum(sum( PolarMatrix(BwN{s}(t-1):BwN{s}(t), BtN(a):BtN(a+1) ) ,2 ));

            L1= sector_sum/area_sector;

            if L1>max_L1
                max_L1=L1;
            end

        end
    end
end

% Tests found that the optimal percentage from the maximum L1 was
% op_thresh_percentage, we find the specific thresh value for this
% image by multiplying this percentage by max_L1
op_thresh_percentage = 0.10;
thresh = max_L1*op_thresh_percentage;

for a=1:(length(BtN)-1)
    for s=2:(length(BwN)-1)
        % Go from the outermost scale towards the center
        for t=length(BwN{s}): -1: 2

            width = BwN{s}(t) -BwN{s}(t-1);
            height = BtN(a+1)-BtN(a);

            area_sector= width*height;

            sector_sum = sum(sum( PolarMatrix(BwN{s}(t-1):BwN{s}(t), BtN(a):BtN(a+1) ) ,2 ));

            L1= sector_sum/area_sector;

            % Making sure all the values set to this one are changed
            % as well
            if L1<thresh
                i=t;
                while i<length(BwN{s})
                    if Scale_lengths{s}(i)==Scale_lengths{s}(i+1)

                        Scale_lengths{s}(i+1)=Bw{s}(t-1);

                    end
                    i=i+1;
                end

                Scale_lengths{s}(t)= Bw{s}(t-1);
            end
        end
    end
end


% Last angle - special case because  the area is larger :
% from last angle detected to the end, plus the area from the
% beginning to the first angle
for t= length(BwN{end}):-1:2

    width = BwN{end}(t) - BwN{end}(t-1);
    height = length(PolarMatrix)- BtN(end) + BtN(1);

    area_sector= width*height;

    sector_sum = sum(sum( PolarMatrix(BwN{end}(t-1):BwN{end}(t), BtN(end):end) ,2 ) ...
        +sum( PolarMatrix(BwN{end}(t-1):BwN{end}(t), 1:BtN(1)) ,2 ) );

    L1= sector_sum/area_sector;

    % Making sure all the values set to this one are changed
    % as well
    if L1<thresh
        i=t;
        while i<length(BwN{end})
            if Scale_lengths{end}(i)==Scale_lengths{end}(i+1)

                Scale_lengths{end}(i+1)=Bw{end}(t-1);

            end
            i=i+1;
        end


        Scale_lengths{end}(t)=Bw{end}(t-1);
    end


end

%Set the old scale values to the new values
Bw = Scale_lengths;

%max_L1

% Fix list of scales because there will be repeated values, this
% loop sets one of the repeated values equal to zero
for s=2:length(Bw)
    for t=1:(length(Bw{s})-1)
        if Bw{s}(t) == Bw{s}(t+1)
            Bw{s}(t)=0;
        end
    end
end


% This for loop gets rid of all the zeros in the list of scales
for s=2:length(Bw)
    Bw{s} = Bw{s}(Bw{s}~=0);
end

% Build the filter bank
mfb = EWT2D_Curvelet_FilterBank(Bw,Bt,W,H,params.option);

% We filter the signal to extract each subband
ff=fft2(f);

ewtc=cell(length(mfb),1);
% We extract the low frequencies first
ewtc{1}=real(ifft2(conj(mfb{1}).*ff));
for s=2:length(mfb)
    ewtc{s}=cell(length(mfb{s}),1);
    for t=1:length(mfb{s})
        if(size(mfb{s}{t},1) > size(ff,1))
            temp = mfb{s}{t};
            mfb{s}{t} = temp((1:size(ff,1)), :);
        end
        if(size(mfb{s}{t},2) > size(ff,2))
            temp = mfb{s}{t};
            mfb{s}{t} = temp(:,(1:size(ff,2)));
        end
        ewtc{s}{t}=real(ifft2(conj(mfb{s}{t}).*ff));
    end
end
