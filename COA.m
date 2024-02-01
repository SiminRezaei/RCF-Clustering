
%%% Designed, Zeinab Montazeri and Developed by Mohammad Dehghani and Pavel Trojovský %%%
function[Best_score,Best_pos,COA_curve]=COA( initial_centers, SearchAgents, Max_iterations, lowerbound, upperbound, dimension, ncluster, na,data)
                                      % COA(lu, nfevalMAX, Np, Nc, initial_centers, data, ncluster, na);
lowerbound=ones(1,dimension).*(lowerbound);                              % Lower limit for variables
upperbound=ones(1,dimension).*(upperbound);                              % Upper limit for variables
%% INITIALIZATION
if initial_centers == 0
    for i=1:dimension
        X(:,i) = lowerbound(i)+rand(SearchAgents,1).*(upperbound(i) - lowerbound(i));                          % Initial population
    end
else
    X = repmat( initial_centers, SearchAgents,1 );
for i =1:SearchAgents
    L=X(i,:);
    fit(i)=FOBJ( data, L, ncluster, na);
end
%%
for t=1:Max_iterations
    %% update the best condidate solution
    [best , location]=min(fit);
    if t==1
        Xbest=X(location,:);                                           % Optimal location
        fbest=best;                                           % The optimization objective function
    elseif best<fbest
        fbest=best;
        Xbest=X(location,:);
    end
    
    
    %%
        for i=1:SearchAgents/2
            
            %% Phase1: Hunting and attacking strategy on iguana (Exploration Phase)
            iguana=Xbest;
            I=round(1+rand(1,1));
            X_P1(i,:)=X(i,:)+rand(1,1) .* (iguana-I.*X(i,:)); % Eq. (4)
            X_P1(i,:) = max(X_P1(i,:),lowerbound);X_P1(i,:) = min(X_P1(i,:),upperbound);
            
            % update position based on Eq (7)
            L=X_P1(i,:);
            F_P1(i)=FOBJ( data, L, ncluster, na);
            if(F_P1(i)<fit(i))
                X(i,:) = X_P1(i,:);
                fit(i) = F_P1(i);
            end
        end
        %%
        for i=1+SearchAgents/2 :SearchAgents
            iguana=lowerbound+rand.*(upperbound-lowerbound); %Eq(5)
            L=iguana;
            F_HL=FOBJ( data, L, ncluster, na);
            I=round(1+rand(1,1));
            
            if fit(i)> F_HL
                X_P1(i,:)=X(i,:)+rand(1,1) .* (iguana-I.*X(i,:)); % Eq. (6)
            else
                X_P1(i,:)=X(i,:)+rand(1,1) .* (X(i,:)-iguana); % Eq. (6)
            end
            X_P1(i,:) = max(X_P1(i,:),lowerbound);X_P1(i,:) = min(X_P1(i,:),upperbound);
            
            % update position based on Eq (7)
            L=X_P1(i,:);
            F_P1(i)=FOBJ( data, L, ncluster, na);
            if(F_P1(i)<fit(i))
                X(i,:) = X_P1(i,:);
                fit(i) = F_P1(i);
            end
        end
        %% END Phase1: Hunting and attacking strategy on iguana (Exploration Phase)
        
        %%
        
        %% Phase2: The process of escaping from predators (Exploitation Phase)
        for i=1:SearchAgents
            LO_LOCAL=lowerbound/t;% Eq(9)
            HI_LOCAL=upperbound/t;% Eq(10)
            
            X_P2(i,:)=X(i,:)+(1-2*rand).* (LO_LOCAL+rand(1,1) .* (HI_LOCAL-LO_LOCAL)); % Eq. (8)
            X_P2(i,:) = max(X_P2(i,:),LO_LOCAL);X_P2(i,:) = min(X_P2(i,:),HI_LOCAL);
            
            % update position based on Eq (11)
            L=X_P2(i,:);
            F_P2(i)=FOBJ( data, L, ncluster, na);
            if(F_P2(i)<fit(i))
                X(i,:) = X_P2(i,:);
                fit(i) = F_P2(i);
            end
            
        end
        %% END Phase2: The process of escaping from predators (Exploitation Phase)
    
    best_so_far(t)=fbest;
    average(t) = mean (fit);
    
end
Best_score=fbest;
Best_pos=Xbest;
COA_curve=best_so_far;
end
