%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function fcenter = mopso(data, np, na, ncluster, tcluster)

%% Problem Definition

% CostFunction=@(x) obj(x);      % Cost Function

nVar = na;             % Number of Decision Variables

VarSize = [ncluster nVar];   % Size of Decision Variables Matrix

VarMin = 0;          % Lower Bound of Variables
VarMax = 1;          % Upper Bound of Variables


%% MOPSO Parameters

MaxIt = 100;           % Maximum Number of Iterations

nPop = 50;           % Population Size

nRep = 50;            % Repository Size

w = 0.5;              % Inertia Weight
wdamp = 0.9;         % Intertia Weight Damping Rate
c1 = 1;               % Personal Learning Coefficient
c2 = 1;               % Global Learning Coefficient

nGrid = 7;            % Number of Grids per Dimension
alpha = 0.1;          % Inflation Rate

beta = 2;             % Leader Selection Pressure
gamma = 2;            % Deletion Selection Pressure

mu = 0.5;             % Mutation Rate

%% Initialization

empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];
empty_particle.IsDominated = [];
empty_particle.GridIndex = [];
empty_particle.GridSubIndex = [];

pop=repmat(empty_particle,nPop,1);

for i=1:nPop
%         sample_id = randi(np, 1, ncluster);
%     pop(i).Position = data(sample_id,:);
    pop(i).Position = unifrnd(VarMin,VarMax,VarSize);

    pop(i).Velocity = zeros(VarSize);
    
    pop(i).Cost = obj(pop(i).Position,data);
    
    
    % Update Personal Best
    pop(i).Best.Position = pop(i).Position;
    pop(i).Best.Cost = pop(i).Cost;
    
end

% Determine Domination
pop = DetermineDomination(pop);

rep = pop(~[pop.IsDominated]);

Grid = CreateGrid(rep,nGrid,alpha);

for i=1:numel(rep)
    rep(i) = FindGridIndex(rep(i),Grid);
end


%% MOPSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        leader=SelectLeader(rep,beta);
        
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
        pop(i).Position = max(pop(i).Position, VarMin);
        pop(i).Position = min(pop(i).Position, VarMax);
        
        pop(i).Cost = obj(pop(i).Position,data);
        
        % Apply Mutation
        pm=(1-(it-1)/(MaxIt-1))^(1/mu);
        if rand<pm
            NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
            NewSol.Cost=obj(NewSol.Position,data);
            if Dominates(NewSol,pop(i))
                pop(i).Position=NewSol.Position;
                pop(i).Cost=NewSol.Cost;

            elseif Dominates(pop(i),NewSol)
                % Do Nothing

            else
                if rand<0.5
                    pop(i).Position=NewSol.Position;
                    pop(i).Cost=NewSol.Cost;
                end
            end
        end
        
        if Dominates(pop(i),pop(i).Best)
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
            
        elseif Dominates(pop(i).Best,pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position=pop(i).Position;
                pop(i).Best.Cost=pop(i).Cost;
            end
        end
        
    end
    
    % Add Non-Dominated Particles to REPOSITORY
    rep=[rep
         pop(~[pop.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep=DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
    
    % Update Grid
    Grid=CreateGrid(rep,nGrid,alpha);

    % Update Grid Indices
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
        
    end
    
    % Plot Costs
%      figure(1);
%      PlotCosts(pop,rep);
%      pause(0.01);
%     
%     % Show Iteration Information
%     disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight
%     wdamp = 1-(it/MaxIt);
     w=w*wdamp;
    
end

%% Select best rep element
rep_size = numel(rep);

for i=1:rep_size
    for j=1:np
        for k=1:ncluster
            diff_c_d = (rep(i).Position(k,:)-data(j,:)).^2;
            dis_p_center1(k) = sqrt(sum(diff_c_d));
        end
        [dis1(j),clu1(j)] = min(dis_p_center1);
    end
    
    for i2=1:ncluster
        ui1 = find(clu1==i2);
        cluster_size1 = size(ui1,2);
        if cluster_size1==0
            [far_point_dis1, far_point1] = max(dis1);
            clu1(far_point1) = i2;
            dis1(far_point1) = pdist2(rep(i).Position(i2,:),data(far_point1,:));  
        end
    end
    
    
    sep_rep(i) = sep(data, ncluster, np, na, clu1, rep(i).Position);
    coh_rep(i) = sum(dis1);
    
    DB_struct = evalclusters(data, clu1', 'DaviesBouldin');
    DB_index(i) = DB_struct.CriterionValues;
    
    sil_dis = silhouette(data, clu1');
    sil_index(i) = sum(sil_dis)./np;
    
    
%     acc1 = confusionmat(tcluster, clu1);
%     acc(i)=sum(diag(acc1))./np;


% [mop, morm mof] = measures(data, np, clu1, tcluster);
% acc(i) = mop;

[mop, morm, mof] = AC(clu1, tcluster', ncluster);
acc(i) = morm;

    
% %     eva = evalclusters(data,clu1','DaviesBouldin');
% %     DB(i) = eva.CriterionValues;
%     
%     s = silhouette(data,clu1);
%     sil_ind(i) = sum(s)./np;
%     
% %     final_ind(i) = ((DB(i).^2)+(sil_ind(i).^2))./2;
%     DB(i) = DBI2(clu1,dis1,rep(i).Position,ncluster);
%     
%     dun(i) = dunns(data, ncluster, clu1, dis1);
end

% min_dun = min(dun);
% max_dun = max(dun);
% dun = (dun-min_dun)./(max_dun-min_dun);

min_sil = min(sil_index);
max_sil = max(sil_index);
sil_index = (sil_index-min_sil)./(max_sil-min_sil);


min_DB = min(DB_index);
max_DB = max(DB_index);
DB_index = (DB_index-min_DB)./(max_DB-min_DB);
DB_index = 1-DB_index;

min_sep = min(sep_rep);
max_sep = max(sep_rep);
sep_rep = (sep_rep-min_sep)./(max_sep-min_sep);

min_coh = min(coh_rep);
max_coh = max(coh_rep);
coh_rep = (coh_rep-min_coh)./(max_coh-min_coh);
coh_rep = 1-coh_rep;

% min_Mi = min(Mi);
% max_Mi = max(Mi);
% Mi = (Mi-min_Mi)./(max_Mi-min_Mi);
% Mi = 1-Mi;

% msil = min(sil_ind);
% max_sil = max(sil_ind);
% sil_ind = (sil_ind-msil)./(max_sil-msil);
% 
% 
% mdun = min(dun);
% max_dun = max(dun);
% Dunn_ind = (dun-mdun)./(max_dun-mdun);
% 
% mdbi = min(DB);
% max_dbi = max(DB);
% DB = (DB-mdbi)./(max_dbi-mdbi);
% DB = 1-DB;
% 
% msil_ind = min(sil_ind);
% max_sil_ind = max(sil_ind);
% sil_ind = (sil_ind-msil_ind)./(max_sil_ind-msil_ind);

final_ind = ((sep_rep.^2)+(coh_rep.^2)+(sil_index).^2)./3;
% final_ind2 = (sep_rep+coh_rep+sil_index)./3;

% [best,ind] = min(DB);

[best, ind] = max(final_ind);
% [best2, ind2] = max(final_ind2);
fcenter = rep(ind).Position;

end
