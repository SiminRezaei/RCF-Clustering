
function z=obj(x,data)
n1=size(x,1);
n2=size(x,2);

np = size(data,1);
nei = round(np*2/100);

for i=1:n1
    for j=1:n1
        diff = (x(i,:)-x(j,:)).^2;
        dist_center(i,j) = sqrt(sum(diff));
    end
    for k=1:np
        diff2 = (x(i,:)-data(k,:)).^2;
        dist2(i,k) = sqrt(sum(diff2));
    end
    sdist2 = sort(dist2(i,:));
    dens(i) = sum(sdist2(1:nei))./nei;
end

% for i=1:n1
%     for j=1:np
%         diff2 = (x(i,:)-data(j,:)).^2;
%         dist2(i,j) = sqrt(sum(diff2));
%     end
% sdist2 = sort(dist2(i,:));
% dens(i) = sum(sdist2(1:10))./9;
% end
% 
% value = min(dist2);
% z3 = sum(value);


ddd = ((n1^2)-n1)./2;
z1 = (sum(sum(dist_center))./2);
 z1 = -z1;
  z2 = sum(dens);
%  z2 = sum(dens);


z = [z1 z2]';

end

