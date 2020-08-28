function strokes_topology = rearrangeLineCoordinates(strokes_topology, vp, vp_indices)
    % Rearranges the coordinates so that the first point is the furtherst from
    % the vanishing point.
    
    mask_orthogonal = (vp_indices < 4);
    
    ind_lines = find(cat(1, strokes_topology(:).primitive_type) == 0);

    ind_orth_lines = ind_lines(mask_orthogonal);
  
    
    lines_coordinates2D = cat(1,strokes_topology(ind_orth_lines).primitive_geom);
    
        
    dist1 = sum((lines_coordinates2D(:, [1,3]) - vp(vp_indices(mask_orthogonal), :)).^2, 2);  
    dist2 = sum((lines_coordinates2D(:, [2,4]) - vp(vp_indices(mask_orthogonal), :)).^2, 2);  

    [~,ind] = max([dist1 dist2],[],2);

    indRearrangeCoordiantes = (ind_orth_lines(ind == 2))';
    
    for i = indRearrangeCoordiantes
        strokes_topology(i).primitive_geom = ...
            strokes_topology(i).primitive_geom([2,1,4,3]);

        strokes_topology(i).points2D(1:end,:) = ...
            strokes_topology(i).points2D(end:-1:1,:);
    end
    
end