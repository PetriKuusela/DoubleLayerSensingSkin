function T = findTriangles(p, g2, H2)

%Use this: https://mathworld.wolfram.com/TriangleInterior.html
v0 = g2(H2(:,1),:);
v1 = g2(H2(:,2),:) - v0;
v2 = g2(H2(:,3),:) - v0;
a = (p(:,1)*v2(:,2)' - p(:,2)*v2(:,1)' - repmat((v0(:,1).*v2(:,2)-v0(:,2).*v2(:,1))',size(p,1),1))./repmat((v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1))',size(p,1),1);
b = -(p(:,1)*v1(:,2)' - p(:,2)*v1(:,1)' - repmat((v0(:,1).*v1(:,2)-v0(:,2).*v1(:,1))',size(p,1),1))./repmat((v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1))',size(p,1),1);
%b = -(p(:,1)*g(H(:,2),2)' - p(:,2)*g(H(:,2),1)' - repmat((g(H(:,1),1).*g(H(:,2),2)-g(H(:,1),2).*g(H(:,2),1))',size(p,1),1))./repmat((g(H(:,2),1).*g(H(:,3),2)-g(H(:,2),2).*g(H(:,3),1))',size(p,1),1);
insides = a>-1e-12 & b>-1e-12 & a+b<1+1e-12;
T = zeros(size(p,1), 1);
for ip = 1:size(p,1)
    temp = find(insides(ip,:));
    if isempty(temp)
        T(ip) = nan;
    else
        T(ip) = temp(1);
    end
end



