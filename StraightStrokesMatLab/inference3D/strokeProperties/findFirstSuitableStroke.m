% Iterates over strokes and searches for the first stroke towards one of
% vanishing points with interction at least with one stroke with dierctions
% towards some other vanishing point.

function si = findFirstSuitableStroke(strokes_topology, intersections)

    si = 1;
    three_dirs = true;
    num_enough_likely = true;
    
    if ismember(strokes_topology(si).line_group, [1,2,3])
        likely_inter_mask = cat(1,intersections(strokes_topology(si).indcs_intrsctns).likely);
        intersecting_strokes = strokes_topology(si).indcs_intrsctng_strks(likely_inter_mask);
        num_enough_likely = length(intersecting_strokes) > 1;
        groups = cat(1,strokes_topology(intersecting_strokes).line_group);
        others = setdiff([1,2,3], strokes_topology(si).line_group);
        three_dirs = ismember(others(1), groups) | ismember(others(2), groups);
    end
    
    while isempty(strokes_topology(si).line_group) || ...
          ~ismember(strokes_topology(si).line_group, [1,2,3])  || ...
          (strokes_topology(si).mean_pressure < 0.05) || ...
          isempty(strokes_topology(si).indcs_intrsctns) || ...
          strokes_topology(si).primitive_type ~= 0  || ~three_dirs || ~num_enough_likely
          
        si = si+1;
        
        if ismember(strokes_topology(si).line_group, [1,2,3])
            likely_inter_mask = cat(1,intersections(strokes_topology(si).indcs_intrsctns).likely);
            intersecting_strokes = strokes_topology(si).indcs_intrsctng_strks(likely_inter_mask);
            num_enough_likely = length(intersecting_strokes) > 1;
            groups = cat(1,strokes_topology(intersecting_strokes).line_group);
            others = setdiff([1,2,3], strokes_topology(si).line_group);
            three_dirs = ismember(others(1), groups) | ismember(others(2), groups);
        end
    end

end
