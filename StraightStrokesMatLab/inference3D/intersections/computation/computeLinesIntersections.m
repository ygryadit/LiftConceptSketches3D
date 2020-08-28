function [points2D, XX, YY] = computeLinesIntersections(lines)
%The points are assumed to lie on the plane that has z = 1
p1 = [lines(:, [1 3]) ones(size(lines, 1), 1)];
p2 = [lines(:, [2 4]) ones(size(lines, 1), 1)];

%Get normals to the planes defined by the line segments 
n = cross(p1, p2, 2);
n = n ./ repmat(sqrt(sum(n.^2,2)), 1, 3);

[XX, YY] = meshgrid(1:size(n,1));

%Compute intersections only once:
A = true(size(n,1));
A = triu(A,1);
XX = XX(A);
YY = YY(A);
ll1 = n(XX(:),:);
ll2 = n(YY(:),:);

%Get the vector that is common for both planes and lie in the image plane
points2D = cross(ll1,ll2);

points2D = [points2D(:,1)./points2D(:,3)...
            points2D(:,2)./points2D(:,3)];
end