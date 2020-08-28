function  [strokes_topology] = computeLength2DLikely(strokes_topology, intersections_to_pair)
for li = 1:length(strokes_topology)
      if strokes_topology(li).primitive_type ~= 0
          strokes_topology(li).length2DLikely = NaN;
          continue;
      end
      strokes_topology(li).length2DLikely = ...
          computeStrokeLengthBetweenIntersectionsApproxLikley(strokes_topology,...
                                                              li,...
                                                              intersections_to_pair);
      
end
end

function length2D = computeStrokeLengthBetweenIntersectionsApproxLikley(strokes_topology, li, intersections)
    if strokes_topology(li).primitive_type ~= 0 || length(strokes_topology(li).indcs_intrsctns) < 2
        length2D = strokes_topology(li).length2DFull;
        return;
    end
    try
    mask_likely = intersections.likely(strokes_topology(li).indcs_intrsctns);
    catch e
        rethrow(e);
    end
    inds_likely = strokes_topology(li).indcs_intrsctns(mask_likely);
    
    nodes__coordinates = intersections.coordinates2D(inds_likely,:);
    
    projected_coordinates = zeros(size(nodes__coordinates,1),1);

    dir = strokes_topology(li).primitive_geom([2,4]) - strokes_topology(li).primitive_geom([1,3]);
    dir = dir./norm(dir);
    for j = 1:size(nodes__coordinates,1)
        projected_coordinates(j,1) = ...
            abs(dot( dir, ...
                     nodes__coordinates(j,:) - ...
                     strokes_topology(li).primitive_geom([1,3]))); 
    
    end
    
    [~, ind_sort] = sort(projected_coordinates);

    if ~isempty(ind_sort)
        length2D = norm(nodes__coordinates(ind_sort(end),:) - nodes__coordinates(ind_sort(1),:));
    else
        length2D = NaN;
    end

end