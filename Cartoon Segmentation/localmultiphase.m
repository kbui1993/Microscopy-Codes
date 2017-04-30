%-------------------------------------------------------------------------
%This function performs local Chan-Vese segmentation onto an image.
%Reference: "An efficient local Chan�Vese model for image segmentation" by
%Wang et. al
%
%Input:
% f - original image
% lambda - weighing parameter for the fidelity term
% beta - weighing parameter for the intensity difference term
% nu - weighing parameter for the regularization term
% tau - time step
%Output:
% u - segmentation result
%-------------------------------------------------------------------------


function [u]=localmultiphase(f, lambda, beta, nu, tau)
% padding parameter
pad=2*lambda;


%creates low-pass filters
[M,N] = size(f);
[XX,YY] = meshgrid((1:N)/N, (1:M)/M);
freqs_1 = XX - 0.5 - 1/N/2;
freqs_2 = YY - 0.5 - 1/M/2;
Lfilter1 = 1 + nu*tau*(freqs_1.^2 + freqs_2.^2);
Lfilter2 = 1 + nu*tau*(freqs_1.^2 + freqs_2.^2);

%creates pre-initialization functions for segmentation (checkerboard functions)
u1 = zeros(M,N);
u2 = zeros(M,N);
for i=1:M
    for j=1:N
        u1(i,j) = sin((pi*i)/3)*sin((pi*j)/3);
        u2(i,j) = sin((pi*i)/10)*sin((pi*j)/10);
    end
end
u1 = double(u1>0);
u2 = double(u2>0);

%padding
f=padarray(f, [pad, pad], 'replicate');
u1=padarray(u1, [pad, pad], 'replicate');
u2=padarray(u2, [pad, pad], 'replicate');
Lfilter1 = padarray(Lfilter1, [pad, pad], 'replicate');
Lfilter2 = padarray(Lfilter2, [pad, pad], 'replicate');

%Filtered image
filtf = imgaussfilt(f,10);

%difference between filtered image and original image
diff = filtf(:) - f(:);
diff1 = filtf - f;

%perform multiphase CV segmentation
for i=1:200

    %Computes mean-value intensities in each regions
    c1=sum(f(:).*(u1(:)==1).*(u2(:)==1))./sum((u1(:)==1).*(u2(:)==1));
    c2=sum(f(:).*(u1(:)==1).*(u2(:)==0))./sum((u1(:)==1).*(u2(:)==0));
    c3=sum(f(:).*(u1(:)==0).*(u2(:)==1))./sum((u1(:)==0).*(u2(:)==1));
    c4=sum(f(:).*(u1(:)==0).*(u2(:)==0))./sum((u1(:)==0).*(u2(:)==0));
    
    d1=sum(diff.*(u1(:)==1).*(u2(:)==1))./sum((u1(:)==1).*(u2(:)==1));
    d2=sum(diff.*(u1(:)==1).*(u2(:)==0))./sum((u1(:)==1).*(u2(:)==0));
    d3=sum(diff.*(u1(:)==0).*(u2(:)==1))./sum((u1(:)==0).*(u2(:)==1));
    d4=sum(diff.*(u1(:)==0).*(u2(:)==0))./sum((u1(:)==0).*(u2(:)==0));
    
    %computing u1
    u1old = u1;
    u1=u1+tau*(-lambda*(c1-f).^2.*u2-lambda*(c2-f).^2.*(1-u2)+lambda*(c3-f).^2.*u2+lambda*(c4-f).^2.*(1-u2));
    u1=u1+tau*(beta*(-(diff1-d1).^2.*u2-(diff1-d2).^2.*(1-u2)+(diff1-d3).^2.*u2+(diff1-d4).^2.*(1-u2)));
    u1=real(ifft2(ifftshift(fftshift(fft2(u1))./Lfilter1)));
    
    %Thresholding u1
    u1 = 1 .* (u1 > 0.5);
    
    %computing u2
    u2old = u2;
    u2=u2+tau*(-lambda*(c1-f).^2.*u1+lambda*(c2-f).^2.*u1-lambda*(c3-f).^2.*(1-u1)+lambda*(c4-f).^2.*(1-u1));
    u2=u2+tau*(beta*(-(diff1-d1).^2.*u1 + (diff1-d2).^2.*u1 - (diff1-d3).^2.*(1-u1)+(diff1-d4).^2.*(1-u1)));

    u2=real(ifft2(ifftshift(fftshift(fft2(u2))./Lfilter2)));

    %Thresholding u2
    u2 = 1.* (u2 > 0.5);
    
    if all(u1old(:) == u1(:)) && all(u2old(:) == u2(:))
        fprintf('Number of iterations completed: %d \n \n', i);
        break;
    end
   
end

if i == 200
    fprintf('Number of iterations completed: %d \n \n', i);
end

%unpadding
u1 = u1(1+pad:end-pad,1+pad:end-pad);
u2 = u2(1+pad:end-pad,1+pad:end-pad);

%combining u1 and u2 to obtain final segmentation result
u1=0.6666*u1;
u2=0.3333*u2;
u=u1+u2;

end

