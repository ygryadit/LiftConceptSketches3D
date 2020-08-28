function inds_new = matchIntersections(intersections_rec, intersections, stroke_topology)
intersections_rec_strs_indcs = cat(1,intersections_rec(:).strokes_indices);
num_inter = length(intersections.coordinates2D) - 1;
inds_new = zeros(num_inter, 1);
for i = 1:num_inter
    strokes_indices = intersections.strokes_indices(i,:);
    strokes_indices = sort(strokes_indices);
    [loca, locb] = ismember(strokes_indices, intersections_rec_strs_indcs, 'rows');
    
    if length(locb) == 1 & loca ~= 0
        inds_new(i) = locb;
    end
    
    if loca == 0
%         ind = min(stroke_topology(strokes_indices(1)).merged_with < strokes_indices(1));
%         if ~isempty(ind)
%             strokes_indices(1) = ind;
%             [loca, locb] = ismember(, intersections_rec_strs_indcs, 'rows');
%             
%         end        
        inds_new(i) = NaN;
    end
    
    if length(locb) > 1
        disp('select closest');
    end
end




end