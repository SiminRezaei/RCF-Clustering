function [dist] = FOBJ( data, centers, ncluster, na)

centers = reshape( centers, ncluster, na);
np = size( data, 1 );

[clu_new, ~, ~] = kmeans(data, ncluster, 'MaxIter', 1, 'start', centers);

dist = 0;
for i = 1:np
    dist = dist + sum((data(i, :) - centers( clu_new(i) )).^2)  ;

end