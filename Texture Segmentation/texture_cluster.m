%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs clustering onto the data matrix obtained from
%obtain_datamatrix.m. First it does k-means using cityblock distance and
%then smooths the result using multiphase MBO method. 
%
%Input:
% data - data matrix
% k - number of clusters
% dt - time step for multiphase MBO
%Output:
% result - clustering result of the image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result, initial] = texture_cluster(data, k, dt)
    %performs k means clustering
    initial = kmeans(data, k, 'dist', 'cityblock', 'Replicates', 10);
    
    %initialization (set 75% of the entries to be zeros)
    initializer = indicator(initial,k);
    ind = randperm(size(initializer,1), ceil(size(initializer,1)*.75));
    initializer(ind,:) = 0;
    
    %performs multiphase MBO
    result = multiclassMBO(data, k, dt, initializer);
    
    %return a hard clustering result
    result = hardcluster(result');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function forms an indciator matrix, where each entry taking value of
%either 0 or 1 indiates whether a data point belongs to a cluster
%corresponding to the index of the column.
%
% Input:
% A - column vector with entries of the cluster number
% Output:
% result - indicator matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = indicator(A, k)

%pre-initialie the indicator matrix
result = zeros(size(A,1), k);

%construct the indicator matrix
for i = 1:k
    result(:,i) = double(A==i);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs a hard clustering on the fuzzy clustering result.
%Input:
% -fuzzydata - membership matrix produced by fuzzy clustering
%Output:
% -result - hard membership matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = hardcluster(fuzzydata)

%finds the maximum probability value and associate the pixel to the cluster
%associated with the maximum probability value
[~,result] = max(fuzzydata, [],1);

%takes transpose of result
result = result';

end