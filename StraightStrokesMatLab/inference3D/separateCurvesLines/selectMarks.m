% function strokes_topology = selectMarks(strokes_topology)
% 
% If the line stroke is identified to be a mark, its primitive type changes
% to 2.

function strokes_topology = selectMarks(strokes_topology)

    h2 = figure(2); 
    hold off;
    set(h2, 'Name', 'Intersection 2D');
    global filepath_sketch_img;
    img = readSketchImg(filepath_sketch_img);
    
    THR_MARK_FRACTION_short = 0.01;
    thr_mark = mean(cat(1, strokes_topology(:).length2DFull))*THR_MARK_FRACTION_short ;
%     vals = strokes_topology(~isnan([strokes_topology(:).accuracy_radius])).accuracy_radius;

%     thr_mark = prctile(vals, 75);
    ind_mark_short = find(cat(1, strokes_topology(:).length2DFull) < thr_mark);
    
    
    
    THR_MARK_FRACTION = 0.2;
    thr_mark = mean(cat(1, strokes_topology(:).length2DFull))*THR_MARK_FRACTION;
    ind_mark = find(cat(1, strokes_topology(:).length2DFull) < thr_mark);
    
    imshow(img);
    hold on;
    for i  = reshape(ind_mark, 1, [])
          plot(cat(1,strokes_topology(i).points2D.x), ...
               cat(1,strokes_topology(i).points2D.y), ...
               'LineWidth',1, 'Color', 'r');
    end
    
    
   
    primitives_types = cat(1, strokes_topology(:).primitive_type);
    ind_lines = find((primitives_types == 0) | (primitives_types == -2));
    
    [mask_mark, ~] = ismember(ind_mark, ind_lines);
    ind_mark = ind_mark(mask_mark);
    ind_mark = ind_mark(cat(1,strokes_topology(ind_mark).line_group) == 4);
    
    ind_mark = unique([ind_mark; ind_mark_short]);
    
    for i  = reshape(ind_mark, 1, [])
          plot(cat(1,strokes_topology(i).points2D.x), ...
               cat(1,strokes_topology(i).points2D.y), ...
               'LineWidth',1, 'Color', 'b');
    end
    hold off;
    
    
    for i = reshape(ind_mark, 1, [])
        strokes_topology(i).primitive_type = 2;
    end
    
    
    
    
    
    
    
end