
function [intersections] = ...
    assignAccuracyRadiusLatestStroke(strokes_topology, intersections)
% Assigns accuracy radius as an accuracy radius of the last stroke.

% Output:    
% 
%   intersections.coordinates2D 
%   intersections.strokes_indices
%   intersections.accuracy_radius
% 
%   points_attraction_to_pair.coordinates2D
%   points_attraction_to_pair.strokes_indices
%   points_attraction_to_pair.accuracy_radius
% 

    % Consider only the intersection with non collinear points:
%     ind_inter_non_collinear = ~intersections(:).collinear;
%     intersections.coordinates2D = intersections.coordinates2D(ind_inter_non_collinear,:);
%     intersections.strokes_indices = intersections.strokes_indices(ind_inter_non_collinear,:);
%     intersections.seg_nums = intersections.seg_nums(ind_inter_non_collinear,:);
%     intersections.p_dist_str_segs = intersections.p_dist_str_segs(ind_inter_non_collinear,:);
%     
%     stroke_indices_inter = intersections(:).strokes_indices;
% %     stroke_indices_inter = max(stroke_indices_inter,[],2);
%     stroke_indices_inter = max(stroke_indices_inter,[],2);
%     
%     intersections.accuracy_radius = cat(1,strokes_topology(stroke_indices_inter).accuracy_radius);
% 
%     % Attraction points
%     ind_attract_non_collinear = find(~points_attraction(:).collinear);
% 
%     points_attraction_to_pair.coordinates2D = points_attraction.coordinates2D(ind_attract_non_collinear,:);
%     points_attraction_to_pair.strokes_indices = points_attraction.strokes_indices(ind_attract_non_collinear,:);
%     
%     points_attraction_to_pair.p_dist_str_segs = points_attraction.p_dist_str_segs(ind_attract_non_collinear,:);
% 
%     stroke_indices_attract = points_attraction_to_pair(:).strokes_indices;
%     stroke_indices_attract = max(stroke_indices_attract,[],2);
% 
%     points_attraction_to_pair.accuracy_radius = cat(1,strokes_topology(stroke_indices_attract).accuracy_radius);

    
    % Accuracy radius:
    stroke_indices_inter = intersections(:).strokes_indices;
    stroke_indices_inter = max(stroke_indices_inter,[],2);    
    intersections.accuracy_radius = cat(1,strokes_topology(stroke_indices_inter).accuracy_radius);

%     % Accuracy radius:
%     stroke_indices_attract = points_attraction_to_pair(:).strokes_indices;
%     stroke_indices_attract = max(stroke_indices_attract,[],2);
%     points_attraction_to_pair.accuracy_radius = cat(1,strokes_topology(stroke_indices_attract).accuracy_radius);
end