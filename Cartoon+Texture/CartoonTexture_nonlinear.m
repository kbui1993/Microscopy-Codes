% -------------------------------------------------------------------------
% CARTOONTEXTURE_NONLINEAR     Cartoon + Texture decomposition (nonlinear)
%
%   [u,v] = CARTOONTEXTURE_NONLINEAR(f, sigma)
%       Decomposes an image into cartoon + texture 
%
%   Input:
%        f        input image
%        sigma    texture scale, mesured in pixel mesh
%
%   Output:
%        u        cartoon image
%        v        texture image
%
% Reference: Buades et al. "Fast cartoon+texture image filters"
% -------------------------------------------------------------------------

function [u,v] = CartoonTexture_nonlinear(f,sigma)

% padding f
[m,n] = size(f);
pad   = 2*sigma;
A = [f(1+pad:-1:2,:); f; f(m-1:-1:m-pad,:)];
f = [A(:,1+pad:-1:2), A, A(:,n-1:-1:n-pad)];

% compute linear cartoon image of f
u = lowPF(f, sigma);

% computing lambda
lambda = (localtv(f, sigma) -localtv(u, sigma))./localtv(f, sigma);

% thresholding lambda
a1 = 0.25; a2 = 0.5;
case1 = lambda<a1;
    lambda(case1) = 0;
case2 = (a1<=lambda) & (lambda<=a2);
    lambda(case2) = (lambda(case2)-a1)/(a2-a1);
case3 = lambda>a2;
    lambda(case3) = 1;

% computing cartoon + texture
u = lambda.*u + (1-lambda).*f;
v = f - u;

% cropping (u,v)
u = u(1+pad:end-pad,1+pad:end-pad);
v = v(1+pad:end-pad,1+pad:end-pad);

end

%% Linear Low Pass Filtering

function [u,v] = lowPF(f,sigma)

% frequency domain
[M,N] = size(f);
[X,Y] = meshgrid((1:N)/N, (1:M)/M);
freqs_1 = X - 0.5 - 1/N/2;
freqs_2 = Y - 0.5 - 1/M/2;

% low-pass filter (L_sigma := 1/filter)
filter = 1 + (2.*pi.*sigma.*(freqs_1.^2 + freqs_2.^2).^(1/2)).^4;

% cartoon image
u = real(ifft2(ifftshift(fftshift(fft2(f))./filter)));

% texture image
v = f - u;
end


%% Local Total Variation

function [u] = localtv(f, sigma)

% compute gradient
[dx, dy] = gradient(f);
gradf = sqrt(dx.^2+dy.^2);

% local total variation
u = lowPF(gradf, sigma);

end
