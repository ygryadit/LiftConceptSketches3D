function  lines_pairs = getPairsCollinearStraightStrokes(intersections)
    
    mask_inter_collinear = logical(intersections.collinear(:));
    lines_pairs = intersections.strokes_indices(mask_inter_collinear,:);
  
    mask_ind = lines_pairs(:,1) > lines_pairs(:,2);
    temp = lines_pairs(mask_ind,1);
    lines_pairs(mask_ind,1) = lines_pairs(mask_ind,2);
    lines_pairs(mask_ind,2) = temp;

    [lines_pairs, ~, ~] = unique(lines_pairs, 'rows');
   
end

