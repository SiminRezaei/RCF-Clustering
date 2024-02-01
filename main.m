
clear;
clc;
close all;
warning('off', 'all');
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


nloop = 30;

%% evaluation of COA algorithm implementation
for i=1:nloop
    disp('-------------- COA kmeans --------------------------')
    fprintf('\n');
    center3 = mopso(data, np, na, ncluster, tcluster);
    initial_centers = reshape( center3,1,  ncluster * na);
    dimension = ncluster * na;
    lowerbound = 0;
    upperbound = 1;
    SearchAgents = 100;
    Max_iterations = 1000;
    
    
%     [ GlobalParams, GlobalMin ] = COA(lu, nfevalMAX, Np, Nc, initial_centers, data, ncluster, na);
    [Best_score,Best_pos,COA_curve]= COA( initial_centers,SearchAgents,Max_iterations,lowerbound,upperbound,dimension, ncluster, na,data);
    %                                COA(lu, nfevalMAX,n_packs,n_coy, initial_centers, data, nclusters, na)
    centers = reshape( Best_pos, ncluster, na );
    [kmeansnovel_clu, kmeansnovel_center, kmeansnovel_MSE] = kmeans(data, ncluster,'MaxIter',1,'Start',centers );
    
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
disp(['COA k means NMI index : ' num2str(kmeansnovel_NMI_mean)]);
disp(['COA k means NMI std : ' num2str(kmeansnovel_NMI_std)]);

fprintf('\n');
disp(['COA k means ARI : ' num2str(kmeansnovel_ARI_mean)]);
disp(['COA k means ARI std : ' num2str(kmeansnovel_ARI_std)]);

fprintf('\n');
disp(['COA k means MS : ' num2str(kmeansnovel_MSscore_mean)]);
disp(['COA k means MS std : ' num2str(kmeansnovel_MSscore_std)]);

%% evaluation of FHO algorithm
% for i=1:nloop
%     disp('-------------- FHO kmeans --------------------------')
%     fprintf('\n');
%     center3 = mopso(data, np, na, ncluster, tcluster);
%     initial_centers = reshape( center3,1,  ncluster * na);
%     VarNumber = ncluster * na;                         % Number of variables;
%     VarMin = zeros(1,VarNumber);            % Lower bound of variable;
%     VarMax = ones(1,VarNumber);             % Upper bound of variable;
%     MaxFes = 1000 ;                       % Maximum number of generations
%     nPop = 25 ;                             % Maximum number of initial candidates
%     HN = randi([1 ceil(nPop/5)],1,1) ;
%     [ Eval_Number, Conv_History, Best_Pos]   = FHO(  VarMin, VarMax, MaxFes, nPop, HN, initial_centers, data, ncluster, na );
%     centers = reshape(Best_Pos, ncluster, na );
%     [kmeansnovel_clu, kmeansnovel_center, kmeansnovel_MSE] = kmeans(data, ncluster,'MaxIter',1,'Start',centers );
%     
%     kmeansnovel_MSE_loop(i) = mse_cal(data, kmeansnovel_center, kmeansnovel_clu, np);
%     [NMI_ind4, accuracy_ind4, precision_ind4] = AC(kmeansnovel_clu', tcluster', ncluster);
%     kmeansnovel_NMI_loop(i) = NMI_ind4;
%     kmeansnovel_ARI_loop(i) = ARI4(kmeansnovel_clu, tcluster);
%     kmeansnovel_MSscore_loop(i) = MSscore(data, ncluster, tcluster, kmeansnovel_clu);
%     fprintf('\n');
%     disp('---------------------------------------------------------------')
% end
% 
% kmeansnovel_NMI_mean = sum(kmeansnovel_NMI_loop)./nloop;
% kmeansnovel_NMI_std = std(kmeansnovel_NMI_loop);
% 
% kmeansnovel_ARI_mean = sum(kmeansnovel_ARI_loop)./nloop;
% kmeansnovel_ARI_std = std(kmeansnovel_ARI_loop);
% 
% kmeansnovel_MSscore_mean = sum(kmeansnovel_MSscore_loop)./nloop;
% kmeansnovel_MSscore_std = std(kmeansnovel_MSscore_loop);
% 
% fprintf('\n');
% disp(['FOH k means NMI index : ' num2str(kmeansnovel_NMI_mean)]);
% disp(['FOH k means NMI std : ' num2str(kmeansnovel_NMI_std)]);
% 
% fprintf('\n');
% disp(['FOH k means ARI : ' num2str(kmeansnovel_ARI_mean)]);
% disp(['FOH k means ARI std : ' num2str(kmeansnovel_ARI_std)]);
% 
% fprintf('\n');
% disp(['FOH k means MS : ' num2str(kmeansnovel_MSscore_mean)]);
% disp(['FOH k means MS std : ' num2str(kmeansnovel_MSscore_std)]);
% 
% %% evaluation of ARO algorithm implementation
% 
%  for i=1:nloop
%     disp('-------------- ARO kmeans --------------------------')
%     fprintf('\n');
%     center3 = mopso(data, np, na, ncluster, tcluster);
%     initial_centers = reshape( center3,1,  ncluster * na);
%     MaxIt = 500;
%     nPop = 25 ; 
%     Dim = ncluster * na; % Maximum number of initial candidates
%     Low = zeros(1, Dim );
%     Up = ones(1, Dim );
%     [ BestX, BestF, HisBestF]  = ARO( Low, Up, Dim, MaxIt, nPop, initial_centers, data, ncluster, na );
%     centers = reshape(BestX, ncluster, na );
%     [kmeansnovel_clu, kmeansnovel_center, kmeansnovel_MSE] = kmeans(data, ncluster,'MaxIter',1,'Start',centers );
%     
%     kmeansnovel_MSE_loop(i) = mse_cal(data, kmeansnovel_center, kmeansnovel_clu, np);
%     [NMI_ind4, accuracy_ind4, precision_ind4] = AC(kmeansnovel_clu', tcluster', ncluster);
%     kmeansnovel_NMI_loop(i) = NMI_ind4;
%     kmeansnovel_ARI_loop(i) = ARI4(kmeansnovel_clu, tcluster);
%     kmeansnovel_MSscore_loop(i) = MSscore(data, ncluster, tcluster, kmeansnovel_clu);
%     fprintf('\n');
%     disp('---------------------------------------------------------------')
% end
% 
% kmeansnovel_NMI_mean = sum(kmeansnovel_NMI_loop)./nloop;
% kmeansnovel_NMI_std = std(kmeansnovel_NMI_loop);
% 
% kmeansnovel_ARI_mean = sum(kmeansnovel_ARI_loop)./nloop;
% kmeansnovel_ARI_std = std(kmeansnovel_ARI_loop);
% 
% kmeansnovel_MSscore_mean = sum(kmeansnovel_MSscore_loop)./nloop;
% kmeansnovel_MSscore_std = std(kmeansnovel_MSscore_loop);
% 
% fprintf('\n');
% disp(['ARO k means NMI index : ' num2str(kmeansnovel_NMI_mean)]);
% disp(['ARO k means NMI std : ' num2str(kmeansnovel_NMI_std)]);
% 
% fprintf('\n');
% disp(['ARO k means ARI : ' num2str(kmeansnovel_ARI_mean)]);
% disp(['ARO k means ARI std : ' num2str(kmeansnovel_ARI_std)]);
% 
% fprintf('\n');
% disp(['ARO k means MS : ' num2str(kmeansnovel_MSscore_mean)]);
% disp(['ARO k means MS std : ' num2str(kmeansnovel_MSscore_std)]);
% 
