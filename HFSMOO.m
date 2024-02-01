function [k_n_clu, k_n_center, k_n_MSE, center3] = HFSMOO(data, ncluster, np, na, tcluster)

center3 = mopso(data, np, na, ncluster, tcluster);

[clu_new, center_new, MSE_novel] = kmeans(data, ncluster,'MaxIter',1,'Start',center3);

k_n_clu = clu_new;

%     for i=1:ncluster
%         ui = find(clu_new==i);
%         cluster_size = size(ui,2);
%         if cluster_size==0
%             [far_point_dis, far_point] = max(dis);
%             clu_new(far_point) = i;
%             ui = far_point;
%             dis1(far_point) = pdist2(center3(i,:),data(far_point,:));
%         end
%     end
    
k_n_center = center_new;
k_n_MSE = MSE_novel;
% iter = 0;
% while iter<1000
%     iter = iter+1;
%     
%     for j=1:np
%         for i=1:ncluster
%             diff_c_d2 = (center3(i,:)-data(j,:)).^2;
%             dis_p_center2(i) = sqrt(sum(diff_c_d2));
%         end
%          [dis(j),clu(j)] = min(dis_p_center2);
%     end
%     
%     for i=1:ncluster
%         ui = find(clu==i);
%         cluster_size = size(ui,2);
%         if cluster_size==0
%             [far_point_dis, far_point] = max(dis);
%             clu(far_point) = i;
%             ui = far_point;
%             dis1(far_point) = pdist2(center3(i,:),data(far_point,:));
%         end
%         up = data(ui,:);
%         if size(ui,2)==1
%             newcenter(i,:) = up;
%         else
%             newcenter(i,:) = mean(up);
%         end
%     end
%     
%     fit(iter) = sum(dis);
%     disp(['Iter =' num2str(iter) '    BEST =' num2str(fit(iter))]);
%     if isequal(newcenter,center3)
%         break
%     else
%         center3 = newcenter;
%     end
% end
% 
% k_n_clu = clu;
% k_n_center = center3;

end

