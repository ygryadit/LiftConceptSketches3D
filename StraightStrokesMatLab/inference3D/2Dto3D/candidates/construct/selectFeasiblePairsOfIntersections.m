function pairs = selectFeasiblePairsOfIntersections(...
                    intersections,...
                    inds_intrsctns_prvs_strks,...
                    pairsInterInter,...
                    cur_stroke_ind,...
                    strokes_topology)
                
    pairs = nchoosek(inds_intrsctns_prvs_strks, 2);              
    pairs = orderIndicesInPairs(pairs);

    % Ignore the pairs that are in 2D are grouping perceptually:
    ind = ismember(pairs, pairsInterInter, 'rows');
    pairs = pairs(~ind, :);


    if isempty(pairs)
        return;
    end

    % Find 2D piece covered by a paper, if less than 10% of the
    % stroke length than do not include such a pair.
%             ind_keep_pairs = find(sqrt(sum(( cat(1,intersections(pairs(:,1)).coordinates2D)- ...
%                                              cat(1,intersections(pairs(:,2)).coordinates2D) ).^2,2))./...
%                                              strokes_topology(cur_stroke.ind).length2DFull > 0.1);
    ind_keep_pairs = find(sqrt(sum(( cat(1,intersections(pairs(:,1)).coordinates2D)- ...
                                     cat(1,intersections(pairs(:,2)).coordinates2D) ).^2,2))./...
                                     strokes_topology(cur_stroke_ind).length2D > 0.1);

    pairs = pairs(ind_keep_pairs,:);

end


function pairs = orderIndicesInPairs(pairs)
    ind = pairs(:,1) > pairs(:,2);
    temp = pairs(ind, 1);
    pairs(ind, 1) = pairs(ind, 2);
    pairs(ind, 2) = temp;
end
