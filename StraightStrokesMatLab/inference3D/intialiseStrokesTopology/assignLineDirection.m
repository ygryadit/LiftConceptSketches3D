function strokes_topology = assignLineDirection(strokes_topology,...
                                           line_group,...
                                           vps_selected,... % additional vanishing points for sets of paralel lines
                                           inds_axis,...    % axis that is orthogonal to the line
                                           inds_active_lines) % indices of lines that belond to a certain cluster)
    primitive_types = cat(1, strokes_topology(:).primitive_type);
%     ind_lines = find((primitive_types == 0) | (primitive_types == -2)| (primitive_types == -3));
    ind_lines = find((primitive_types == 0));
    
    for i = 1:length(ind_lines)
        
        ind_line = ind_lines(i);        
        strokes_topology(ind_line).line_group = line_group(i);
        
    end
    
    for i = 1:length(inds_active_lines)
        for ii = 1:length(inds_active_lines{i})
            strokes_topology(ind_lines(inds_active_lines{i}(ii))).line_group = 5;
            strokes_topology(ind_lines(inds_active_lines{i}(ii))).ind_orth_ax = inds_axis(i);
        end
    end
    
    
end
