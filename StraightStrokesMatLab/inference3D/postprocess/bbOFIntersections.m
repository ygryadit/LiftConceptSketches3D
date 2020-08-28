function [bound_box, intersections_active_coord3D, inds] = bbOFIntersections(intersections)

for i = 1:length(intersections)
    if isempty(intersections(i).is_active)
        intersections(i).is_active = 0;
    end
end

inds = find(( [intersections(:).is_active] == 1));

intersections_active = intersections(cat(1,intersections(:).is_active) == 1);
mask = cat(1,intersections_active(:).collinear) == 0;

inds  = inds(mask);


% intersections_active = intersections(cat(1,intersections(:).is_active) == 1);
% intersections_active = intersections_active(cat(1,intersections_active(:).collinear) == 0); 
intersections_active_coord3D = cat(1,intersections(inds).coordinates3D);
if isempty(intersections_active_coord3D)
    return;
end
bound_box(1,1) = min(intersections_active_coord3D(:,1));
bound_box(1,2) = max(intersections_active_coord3D(:,1));
bound_box(2,1) = min(intersections_active_coord3D(:,2));
bound_box(2,2) = max(intersections_active_coord3D(:,2));
bound_box(3,1) = min(intersections_active_coord3D(:,3));
bound_box(3,2) = max(intersections_active_coord3D(:,3));

p = 0.005;
dx = p*abs(bound_box(1,2) - bound_box(1,1));
dy = p*abs(bound_box(2,2) - bound_box(2,1));
dz = p*abs(bound_box(3,2) - bound_box(3,1));


bound_box(1,1) = bound_box(1,1) - dx;
bound_box(1,2) = bound_box(1,2) + dx;

bound_box(2,1) = bound_box(2,1) - dy;
bound_box(2,2) = bound_box(2,2) + dy;

bound_box(3,1) = bound_box(3,1) - dz;
bound_box(3,2) = bound_box(3,2) + dz;

end