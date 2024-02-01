function mse_val = mse_cal(data, center_ind, clu_ind, np)

for i=1:np
    dist_center(i) = (pdist2(data(i,:),center_ind(clu_ind(i),:))).^2;
end
mse_val = sum(dist_center);
end
