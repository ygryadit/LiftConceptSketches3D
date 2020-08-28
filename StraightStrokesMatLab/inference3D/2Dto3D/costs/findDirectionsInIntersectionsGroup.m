function dirs = findDirectionsInIntersectionsGroup(inds_cls_intrsctns,...
                                            inds_intrsctns__assigned,...
                                            inds_intrsctns__mult_cnddts,...
                                            inds_intrsctns__mult_cnddts_ind,...
                                            intersections,...
                                            strokes_topology,...
                                            cur_stroke_ind)

    dirs = zeros(length(inds_cls_intrsctns), 3);   
    %% Assigned
    inds_cls_intrsctns_assigned = intersect(inds_cls_intrsctns, inds_intrsctns__assigned);

    strokes_at_inter = cat(1,intersections(inds_cls_intrsctns_assigned).strokes_indices);
    strokes_at_inter = strokes_at_inter(strokes_at_inter ~= cur_stroke_ind);

    num_assigned = length(inds_cls_intrsctns_assigned);
    for i= 1:num_assigned 
        dirs(i,:) = strokes_topology(strokes_at_inter(i)).direction_vec;
    end

    %% Multiple
    [inds_cls_intrsctns_mult, ~, locs] = intersect(inds_cls_intrsctns, inds_intrsctns__mult_cnddts);
    if isempty(inds_cls_intrsctns_mult)
        return;
    end
    
    inds_cls_intrsctns_mult_ind = inds_intrsctns__mult_cnddts_ind(locs);

%     strokes_at_inter2 = cat(1,intersections(inds_cls_intrsctns_mult).strokes_indices);
    
%     strokes_at_inter2 =  strokes_at_inter2(strokes_at_inter2 ~= cur_stroke_ind);

    for i= 1:length(inds_cls_intrsctns_mult)    
        mask_pair = intersections(inds_cls_intrsctns_mult(i)).strokes_indices ~= cur_stroke_ind;
        stroke_at_inter2 = intersections(inds_cls_intrsctns_mult(i)).strokes_indices(mask_pair);
        inds_cnddt_ln = intersections(inds_cls_intrsctns_mult(i)).cnddts3D(inds_cls_intrsctns_mult_ind(i)).cnddt_lns{mask_pair};
        dirs(i+num_assigned,:) = strokes_topology(stroke_at_inter2).candidate_lines(inds_cnddt_ln).dir;
    end

end