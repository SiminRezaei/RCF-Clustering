
clear all;
clc;
close all;

% % 
% distcomp.feature('LocalUseMpiexec',false)
% poolobj=gcp('nocreate');
% if ~isempty(poolobj)
%     delete(poolobj);
% end
% parpool('local',24);


load('iris.mat');
data = iris(:,1:end-1);
tcluster = iris(:,end);
ncluster = max(tcluster);
np = size(data,1);
na = size(data,2);

for i=1:na
    mina1 = min(data(:,i));
    maxa1 = max(data(:,i));
    if mina1==maxa1
        continue;
    else
        data(:,i) = (data(:,i)-mina1)/(maxa1-mina1);
    end
end


nloop = 1;

for i=1:nloop
    disp('-------------- Novel kmeans Itration --------------------------')
    fprintf('\n');
    [kmeansnovel_clu, kmeansnovel_center, kmeansnovel_MSE, centers3] = HFSMOO(data, ncluster, np, na, tcluster);
    kmeansnovel_MSE_loop(i) = mse_cal(data, kmeansnovel_center, kmeansnovel_clu, np);
    [NMI_ind4, accuracy_ind4, precision_ind4] = AC(kmeansnovel_clu', tcluster', ncluster);
    kmeansnovel_NMI_loop(i) = NMI_ind4;
    kmeansnovel_ARI_loop(i) = ARI4(kmeansnovel_clu, tcluster);
    kmeansnovel_MSscore_loop(i) = MSscore(data, ncluster, tcluster, kmeansnovel_clu);
    fprintf('\n');
    disp('---------------------------------------------------------------')
end

kmeansnovel_NMI_mean = sum(kmeansnovel_NMI_loop)./nloop;
kmeansnovel_NMI_std = std(kmeansnovel_NMI_loop);

kmeansnovel_ARI_mean = sum(kmeansnovel_ARI_loop)./nloop;
kmeansnovel_ARI_std = std(kmeansnovel_ARI_loop);

kmeansnovel_MSscore_mean = sum(kmeansnovel_MSscore_loop)./nloop;
kmeansnovel_MSscore_std = std(kmeansnovel_MSscore_loop);

fprintf('\n');
disp(['Novel k means NMI index : ' num2str(kmeansnovel_NMI_mean)]);
disp(['Novel k means NMI std : ' num2str(kmeansnovel_NMI_std)]);

fprintf('\n');
disp(['Novel k means ARI : ' num2str(kmeansnovel_ARI_mean)]);
disp(['Novel k means ARI std : ' num2str(kmeansnovel_ARI_std)]);

fprintf('\n');
disp(['Novel k means MS : ' num2str(kmeansnovel_MSscore_mean)]);
disp(['Novel k means MS std : ' num2str(kmeansnovel_MSscore_std)]);



