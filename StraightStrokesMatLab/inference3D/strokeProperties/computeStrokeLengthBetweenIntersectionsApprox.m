function length2D = computeStrokeLengthBetweenIntersectionsApprox(strokes_topology, li, intersections)
    if strokes_topology(li).primitive_type ~= 0 || length(strokes_topology(li).indcs_intrsctns) < 2
        length2D = strokes_topology(li).length2DFull;
        return;
    end
    
    nodes__coordinates = intersections.coordinates2D((strokes_topology(li).indcs_intrsctns),:);
    
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
    
    length2D = norm(nodes__coordinates(ind_sort(end),:) - nodes__coordinates(ind_sort(1),:));
    
end