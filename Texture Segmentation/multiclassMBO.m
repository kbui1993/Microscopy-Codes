%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MULTIPHASEMBO     MBO clustering.
%
%   MULTIPHASEMBO(D,K,dt,U_hat)
%       This function performs MBO clustering on a data image.
%
%   Input:
%        -D       data matrix
%        -K       number of clusters
%        -dt      time step
%        -U_hat   initialization
%
%   Output:
%        -U       clustering result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = multiclassMBO(D, K, dt, U_hat)

%number of times to threshold
Ns = 3;

%weighing parameter for fidelity term
mu = 30;

%obtaining the eigenvectors and eigenvalues using Nystrom Extension
[X, Lambda] = nystrom_colo2(D, 300, 30);

%obtaining size of matrix
[~, n] = size(X); 

%calculating important variable
Y = (eye(n) +(dt/Ns)*Lambda)\transpose(X);

%initialization step for U
U = U_hat;
Uold = U_hat;

%tolerance parameter
nu = 10^(-7);
tol = 10000;

%compute cluster
while tol > nu
    for k = 1:K
    
        %MBO Step
        for i = 1:Ns
            Z = Y*(U-(dt/Ns)*mu*(U-U_hat));
            U = double(X*Z);
        end
        
        I = eye(K);
    
    
        %Thresholding by projecting to simplex
        U = double(projsplx(U));
        dist = pdist2(U,I, 'cityblock');
        [~,pos] = min(dist,[],2);
        U = I(pos,:);
    end
    %compute tolerance
    diffU = U-Uold;
    diffUnorm = diag(sqrt(diffU'*diffU));
    Uoldnorm = diag(sqrt(Uold'*Uold));
    tol = max(diffUnorm)/max(Uoldnorm);
    
    %update Uold
    Uold = U;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project an n-dim vector y to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
%
%The original algorithm is provided by
% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
%
% Jan. 14, 2011.
%
%It has been modified so that it runs more quickly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = projsplx(y)


m = size(y,2); 
n = size(y,1);
s = sort(y, 2, 'descend'); tmpsum = zeros(n,1);

ind = zeros(n,1);
new_ind = zeros(n,1);
tmax = zeros(n,1);
for ii = 1:m-1
    tmpsum = tmpsum + s(:,ii).*(1-ind);
    tmax = tmax.*ind+(1-ind).*(tmpsum - 1)/ii;

    new_ind = new_ind + (tmax>=s(:,ii+1));
    ind = new_ind>0;
    
end

tmax2 = (tmpsum + s(:,m) -1)/m;

tmax3 = ind.*tmax + tmax2.*(1-ind);

x = max(bsxfun(@minus, y, tmax3), 0);

return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by REU 2012 Group
%This functionn performs Nystrom Extension version of spectral clustering.
% Input: 
% data: In vector form
% num_samples: number of eigenvalues/vectors, also the number of random
% samples
% sigma: Adjusts the gaussian similarity
%
% Output: 
% V: eigenvectors
% L: eigenvalues of Laplacian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V, L] = nystrom_colo2(data, num_samples, sigma) 

% randomly select samples
num_rows = size(data, 1);
permed_index = randperm(num_rows);
sample_data = data(permed_index(1:num_samples), :, :);
other_data = data(permed_index(num_samples+1:num_rows), :, :);
clear data;

% calculate the distance between samples themselves
A = pdist2(sample_data, sample_data, 'cityblock');

A = exp(-A/sigma);


% calculate the distance between samples and other points
other_points = num_rows - num_samples;
B = pdist2(sample_data, other_data, 'cityblock');

B = exp(-B/sigma);

% clear sample_data other_data;


% Normalize A and B using row sums of W, where W = [A B; B' B'*A^-1*B].
% Let d1 = [A B]*1, d2 = [B' B'*A^-1*B]*1, dhat = sqrt(1./[d1; d2]).
B_T = B';
d1 = sum(A, 2) + sum(B, 2);
d2 = sum(B_T, 2) + B_T*(pinv(A)*sum(B, 2));
dhat = sqrt(1./[d1; d2]);
A = A .* (dhat(1:num_samples)*dhat(1:num_samples)');
B1 = dhat(1:num_samples)*dhat(num_samples+(1:other_points))';
B = B .* B1;
% clear d1 d2 B1 dhat;

% Do orthogalization and eigendecomposition
Asi = sqrtm(pinv(A));
B_T = B';
BBT = B*B_T;
W = single(zeros(size(A, 1)+size(B_T, 1), size(A, 2)));
W(1:size(A, 1), :) = A;
W(size(A, 1)+1:size(W, 1), :) = B_T;
% clear B B_T;
% Calculate R = A + A^-1/2*B*B'*A^-1/2
R = A + Asi*BBT*Asi;
R = (R + R')/2; % Make sure R is symmetric, sometimes R can be non-symmetric because of numerical inaccuracy
[U, L] = eigs(R,size(R,2));
[~, ind] = sort(diag(L), 'descend');
U = U(:, ind); % in decreasing order
L = L(ind, ind); % in decreasing order
clear A R BBT;
W = W*Asi;
V = W*U(:, 1:num_samples)*pinv(sqrt(L(1:num_samples, 1:num_samples)));

V(permed_index,:) = V;
V = real(V);
L = 1-diag(L);
L = diag(L);
end
