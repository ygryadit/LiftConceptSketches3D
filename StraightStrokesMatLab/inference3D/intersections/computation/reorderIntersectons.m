% Orders intersection in the order of their appearane in the sketch
% Keep only likely ones
function intersections = reorderIntersectons(intersections)
    [intersections.strokes_indices, ind]    = sortrows(intersections.strokes_indices);
    vals = intersections.strokes_indices(:,1) + intersections.strokes_indices(:,2);
    [~, ind2] = sort(vals);
    ind = ind(ind2);
    intersections.strokes_indices = intersections.strokes_indices(ind2,:);
    
    
    
    intersections.coordinates2D             = intersections.coordinates2D(ind, :);
    intersections.collinear                 = intersections.collinear(ind, :);
    intersections.seg_nums                  = intersections.seg_nums(ind, :);
    intersections.p_dist_str_segs           = intersections.p_dist_str_segs(ind, :);
    intersections.accuracy_radius           = intersections.accuracy_radius(ind, :);
   
    if isfield(intersections, 'strokes_indices_original')
        intersections.strokes_indices_original  = intersections.strokes_indices_original(ind, :);
    end
    if isfield(intersections, 'likely')
        intersections.likely = intersections.likely(ind);
    end
    if isfield(intersections, 'tangent')
        intersections.tangent = intersections.tangent(ind);
    end
    
%     intersections.likely                    = intersections.likely(ind);
%     
%     % keep only likely ones:
%     ind_likely = find(intersections.likely);
%     intersections.strokes_indices = intersections.strokes_indices(ind_likely,:);
%     intersections.coordinates2D             = intersections.coordinates2D(ind_likely, :);
%     intersections.p_dist_str_segs           = intersections.p_dist_str_segs(ind_likely, :);
%     intersections.accuracy_radius           = intersections.accuracy_radius(ind_likely, :);
%     intersections.likely                    = intersections.likely(ind_likely);


end