function seperation = sep(data, ncluster, np, na, clu, center_i)

data_center = (sum(data)./np);

% for i=1:ncluster
%     clu_i = find(clu==i);
%     n_i(i) = size(clu_i,1);
%     data_i = data(clu_i,:);
%     s_data_i = size(data_i,1);
%     if s_data_i==1
%         center_i(i,:) = data(clu_i,:);
%     else
%         center_i(i,:) = (sum(data_i)./n_i(i));
%     end
% end

for i=1:ncluster
    clu_i = find(clu==i);
    n_i(i) = size(clu_i,2);
end

for i=1:ncluster
    dist_center(i) = sqrt(sum((center_i(i,:)-data_center).^2));
    dist_center(i) = dist_center(i)*n_i(i);
end

seperation = sum(dist_center);
end