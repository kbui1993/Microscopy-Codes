%------------------------------------------------------------------------
%This function performs multiphase CV segmentation using the MBO Scheme. In
%other words, this segments an image up to 4 regions. The energy functional
%is the following:
%
% E_{CV} = \lambda_1 \int_{\Omega} (\mu_1 - I_0)^2 H(\phi_1)H(\phi_2) +
%          \lambda_2 \int_{\Omega} (\mu_2 - I_0)^2 H(\phi_1)(1-H(\phi_2)) + 
%          \lambda_3 \int_{\Omega} (\mu_3 - I_0)^2 (1-H(\phi_1))H(\phi_2) +
%          \lambda_4 \int_{\Omega} (\mu_4 - I_0)^2 (1-H(\phi_1))(1-H(\phi_2)) + 
%          \beta_1 \int_{\Omega} | \nabla H(\phi_1)| + 
%          \beta_2 \int_{\Omega} | \nabla H(\phi_2)|
%
%Input:
% -f - input image
% -lambda - weighing parameter for fidelity term
% -mu - weighing parameter for regularization term 
% -dt - time step
%Output:
% -u - segmentation result
%--------------------------------------------------------------------------
function [u]=multiphase(f, lambda, mu, dt)

%padding parameter
pad=2*lambda;

% obtains image size
[M,N] = size(f);

%creates low-pass filters
[XX,YY] = meshgrid((1:N)/N, (1:M)/M);
freqs_1 = XX - 0.5 - 1/N;
freqs_2 = YY - 0.5 - 1/M;
Lfilter1 = 1 + mu*dt*(freqs_1.^2 + freqs_2.^2);
Lfilter2 = 1 + mu*dt*(freqs_1.^2 + freqs_2.^2);

%creates pre-initialization functions for segmentation (square functions)
u1 = zeros(M,N);
u2 = zeros(M,N);
for i=1:M
    for j=1:N
        u1(i,j) = sin((pi*i)/3)*sin((pi*j)/3);
        u2(i,j) = sin((pi*i)/10)*sin((pi*j)/10);
    end
end

%threhold
u1 = double(u1>0);
u2 = double(u2>0);

%padding
f=padarray(f, [pad, pad], 'replicate');
u1=padarray(u1, [pad, pad], 'replicate');
u2=padarray(u2, [pad, pad], 'replicate');
Lfilter1 = padarray(Lfilter1, [pad, pad], 'replicate');
Lfilter2 = padarray(Lfilter2, [pad, pad], 'replicate');

%create a matrix of small values
[m,n]= size(u1);
error = 10^(-5)*eye(m,n);

%perform multiphase CV segmentation
for i=1:200

    %Computes mean-value intensities in each regions
    c1=sum(f(:).*(u1(:)==1).*(u2(:)==1))./sum((u1(:)==1).*(u2(:)==1)+error(:));
    c2=sum(f(:).*(u1(:)==1).*(u2(:)==0))./sum((u1(:)==1).*(u2(:)==0)+error(:));
    c3=sum(f(:).*(u1(:)==0).*(u2(:)==1))./sum((u1(:)==0).*(u2(:)==1)+error(:));
    c4=sum(f(:).*(u1(:)==0).*(u2(:)==0))./sum((u1(:)==0).*(u2(:)==0)+error(:));
    
    %computing u1
    u1old = u1;
    u1=u1+dt*(-lambda*(c1-f).^2.*u2-lambda*(c2-f).^2.*(1-u2)+lambda*(c3-f).^2.*u2+lambda*(c4-f).^2.*(1-u2));
    u1=real(ifft2(ifftshift(fftshift(fft2(u1))./Lfilter1)));
    
    %Thresholding u1
    u1 = 1 .* (u1 > 0.50);
    
    %computing u2
    u2old = u2;
    u2=u2+dt*(-lambda*(c1-f).^2.*u1+lambda*(c2-f).^2.*u1-lambda*(c3-f).^2.*(1-u1)+lambda*(c4-f).^2.*(1-u1));
    u2=real(ifft2(ifftshift(fftshift(fft2(u2))./Lfilter2)));

    %Thresholding u2
    u2 = 1.* (u2 > 0.50);
    
    %stopping criterion
    if all(u1old(:) == u1(:)) && all(u2old(:) == u2(:))
        fprintf('Number of iterations completed: %d \n \n', i);
        break;
    end
   
end

%print if the algorithm does not converge
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
