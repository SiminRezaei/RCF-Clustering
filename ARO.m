
function [BestX,BestF,HisBestF]=ARO( Low ,Up ,Dim ,MaxIt,nPop, initial_centers, data, ncluster, na )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FunIndex: Index of function.                       %
    % MaxIt: Maximum number of iterations.               %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PopPos=zeros(nPop,Dim);
PopFit=zeros(nPop,1);
for i=1:nPop
    PopPos(i,:)= initial_centers;
    PopFit(i)= FOBJ(data, PopPos(i,:), ncluster, na );
end
BestF=inf;
BestX=[];
for i=1:nPop
    if PopFit(i)<=BestF
        BestF=PopFit(i);
        BestX=PopPos(i,:);
    end
end
HisBestF=zeros(MaxIt,1);
for It=1:MaxIt
    Direct1=zeros(nPop,Dim);
    Direct2=zeros(nPop,Dim);
    theta=2*(1-It/MaxIt);
    for i=1:nPop
        L=(exp(1)-exp(((It-1)/MaxIt)^2))*(sin(2*pi*rand));%Eq.(3)
        rd=ceil(rand*(Dim));
        Direct1(i,randperm(Dim,rd))=1;
        c=Direct1(i,:);%Eq.(4)
        R=L.*c;%Eq.(2)
        
        A=2*log(1/rand)*theta;% Eq.(15)
        if A>1
            K=[1:i-1 i+1:nPop];
            RandInd=K(randi([1 nPop-1]));
            newPopPos=PopPos(RandInd,:)+R.*( PopPos(i,:)-PopPos(RandInd,:))...
                +round(0.5*(0.05+rand))*randn;%Eq.(1)
        else
            Direct2(i,ceil(rand*Dim))=1;
            gr=Direct2(i,:);%Eq.(12)
            H=((MaxIt-It+1)/MaxIt)*randn; %Eq.(8)
            b=PopPos(i,:)+H*gr.*PopPos(i,:);%Eq.(13)
            newPopPos=PopPos(i,:)+ R.*(rand*b-PopPos(i,:));%Eq.(11)
        end
        newPopPos=SpaceBound(newPopPos,Up,Low);
        newPopFit=FOBJ(data, newPopPos, ncluster, na );
        if newPopFit<PopFit(i)
            PopFit(i)=newPopFit;
            PopPos(i,:)=newPopPos;
        end
    end
    for i=1:nPop
        if PopFit(i)<BestF
            BestF=PopFit(i);
            BestX=PopPos(i,:);
        end
    end
    HisBestF(It)=BestF;
end
