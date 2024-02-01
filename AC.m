function [NMIind, accuracyind, precisionind] = AC(u, v, ncluster)

% u: results of clustering algorithm
% v: true clustering
%
% u = [1 1 1 1 1 2 3 3 1 2 2 2 2 2 3 3 3];
% v = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3];
% ncluster = 3;

np = size(u, 2);
ku = max(u);
kv = max(v);
m = zeros(ku, kv);
for i=1:np
    m(u(i), v(i)) = m(u(i), v(i))+1;
end
mu = sum(m, 2);
mv = sum(m, 1);

vector = max(m');
vectorp = vector./mu';

for i=1:ncluster
    for j=1:ncluster
        if m(i,j)==0
            newmat(i,j) = 0;
        else
            newmat(i,j) = ((m(i,j)./np)*(log2((m(i,j)*np)./(mu(i)*mv(j)))));
            
        end
    end
end

up = 2*sum(sum(newmat));
mu_size = size(mu,1);
for i=1:mu_size
    if mu(i)==0
        newmu(i)=0;
    else
        newmu(i) = (mu(i)./np).*(log2(mu(i)./np));
    end
end
newmu = -newmu;
v1 = sum(newmu);
mv_size = size(mv,2);
for i=1:mv_size
    if mv(i)==0
        newmv(i) = 0;
    else
        newmv(i) = (mv(i)./np).*(log2(mv(i)./np));
    end
end
newmv = -newmv;
v2 = sum(newmv);
NMIind = up./(v1+v2);
accuracyind = (sum(vector)./np);
precisionind = sum(vectorp)./ncluster;

end




