
function [MS] = MSscore(data,ncluster,tcluster,clu)

np = size(data,1);
na = size(data,2);

for i=1:np
    for j=1:np
        if clu(i)==clu(j)
            M(i,j) = 1;
        else
            M(i,j) = 0;
        end
    end
end

for i=1:np
    for j=1:np
        if tcluster(i)==tcluster(j)
            T(i,j) = 1;
        else
            T(i,j) = 0;
        end
    end
end

TM = T-M;

TMnorm = sqrt(sum(sum(TM.^2)));
Tnorm = sqrt(sum(sum((T.^2))));

MS = TMnorm/Tnorm;
end
