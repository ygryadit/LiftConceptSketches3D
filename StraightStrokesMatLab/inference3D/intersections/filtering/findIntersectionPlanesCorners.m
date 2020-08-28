% function findIntersectionPlanesCorners()
% 
% Description
% Return indices of intersections that are corners of planes in 2d: corners
% of cycles of degree 4, there opposite strokes are paralle line.

function inds_corners = findIntersectionPlanesCorners(intersections, strokes_topology)
inds_strks_1 = [intersections.strokes_indices(:,1)];
inds_strks_2 = [intersections.strokes_indices(:,2)];
% Select intersections that are between lines towards vp:

inds_vp12 = find(([strokes_topology(inds_strks_1).line_group] == 1) & ...
     ([strokes_topology(inds_strks_2).line_group] == 2));
inds_vp21 = find(([strokes_topology(inds_strks_1).line_group] == 2) & ...
     ([strokes_topology(inds_strks_2).line_group] == 1)); 

inds_vp13 = find(([strokes_topology(inds_strks_1).line_group] == 1) & ...
     ([strokes_topology(inds_strks_2).line_group] == 3));
inds_vp31 = find(([strokes_topology(inds_strks_1).line_group] == 3) & ...
     ([strokes_topology(inds_strks_2).line_group] == 1));

inds_vp23 = find(([strokes_topology(inds_strks_1).line_group] == 2) & ...
     ([strokes_topology(inds_strks_2).line_group] == 3));
inds_vp32 = find(([strokes_topology(inds_strks_1).line_group] == 3) & ...
     ([strokes_topology(inds_strks_2).line_group] == 2));

inds_corners = [];

inds_corners = getIndsCorners(inds_vp12, ...
                            inds_corners, ...
                            1, ...
                            2, ...
                            inds_strks_1,...
                            inds_strks_2,...
                            strokes_topology,...
                            intersections);
                        
inds_corners = getIndsCorners(inds_vp21, inds_corners, 2, 1, inds_strks_1, inds_strks_2,strokes_topology, intersections);
inds_corners = getIndsCorners(inds_vp13, inds_corners, 1, 3, inds_strks_1, inds_strks_2,strokes_topology, intersections);
inds_corners = getIndsCorners(inds_vp31, inds_corners, 3, 1, inds_strks_1, inds_strks_2,strokes_topology, intersections);
inds_corners = getIndsCorners(inds_vp23, inds_corners, 2, 3, inds_strks_1, inds_strks_2,strokes_topology, intersections);
inds_corners = getIndsCorners(inds_vp32, inds_corners, 3, 2, inds_strks_1, inds_strks_2,strokes_topology, intersections);

end


