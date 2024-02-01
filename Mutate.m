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

function xnew=Mutate(x,pm,VarMin,VarMax)

[row nVar]=size(x);
% for i=1:row
i=randi([1 row]);
j=randi([1 nVar]);

dx=pm*(VarMax-VarMin);

lb=x(i,j)-dx;
if lb<VarMin
    lb=VarMin;
end

ub=x(i,j)+dx;
if ub>VarMax
    ub=VarMax;
end

xnew=x;
xnew(i,j)=unifrnd(lb,ub);
% end

end