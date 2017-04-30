%-------------------------------------------------------------------------
%This function performs local Chan-Vese segmentation onto an image.
%Reference: "An efficient local Chan�Vese model for image segmentation" by
%Wang et. al
%
%Input:
% f - original image
% lambda - weighing parameter for the fidelity term
% beta - weighing parameter for the regularization term
% nu - weighing parameter for the intensity difference term
% tau - time step
%Output:
% u - segmentation result
%-------------------------------------------------------------------------

function [u] = twophase_localmultiphase(f, lambda, beta, nu, tau)

%padding parameter
pad = 2*lambda;

%creates low-pass filters
[M,N] = size(f);
[XX,YY] = meshgrid((1:N)/N, (1:M)/M);
freqs_1 = XX - 0.5 - 1/N/2;
freqs_2 = YY - 0.5 - 1/M/2;
Lfilter = 1 + nu*tau*(freqs_1.^2 + freqs_2.^2);

%creates pre-initialization function for segmentation
u = zeros(M,N);
for i=1:M
    for j=1:N
        u(i,j) = sin((pi*i)/3)*sin((pi*j)/3);
    end
end
u = double(u>0);

%padding
f=padarray(f, [pad, pad], 'replicate');
u=padarray(u, [pad, pad], 'replicate');
Lfilter = padarray(Lfilter, [pad, pad], 'replicate');

%Filtered image
filtf = imgaussfilt(f,10);

%difference between filtered image and original image
diff = filtf(:) - f(:);
diff1 = filtf - f;

%perform CV segmentation
for i = 1:200
    %compute mean value intensities in each region
    c1 = sum(f(:).*(u(:)==1))./sum(u(:)==1);
    c2 = sum(f(:).*(u(:)==0))./sum(u(:)==0);
    
    %computing mean difference between filtered image and original image
    d1 = sum(diff.*(u(:)==1))./sum(u(:)==1);
    d2 = sum(diff.*(u(:)==0))./sum(u(:)==0);
    
    %updating u
    uold = u;
    u = u-tau*lambda*((c1-f).^2-(c2-f).^2);
    u = u -tau*beta*((diff1-d1).^2-(diff1-d2).^2);
    
    u=real(ifft2(ifftshift(fftshift(fft2(u))./Lfilter)));
    
    %Thresholding u
    u = 1.*(u > 0.5);
    
    %stopping criterion
    if all(uold(:) == u(:))
        fprintf('Number of iterations completed: %d \n \n', i);
        break;
    end
    
end

if i == 200
    fprintf('Number of iterations completed: %d \n \n', i);
end


%unpadding
u = u(1+pad:end-pad,1+pad:end-pad);

end