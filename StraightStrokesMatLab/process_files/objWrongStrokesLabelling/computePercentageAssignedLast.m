function percentage_assigned_last = computePercentageAssignedLast(set_removed, strokes_topology_rough)
    
    strokes_analyse = strokes_topology_rough(set_removed);
   
    last_stroke = max([strokes_topology_rough(:).created]);
    
    
    indices_known_creation = [];
    for i = 1:length(strokes_analyse)
        if ~isempty(strokes_analyse(i).created)
            indices_known_creation(end+1) = i;
        end
    end
    
    strokes_analyse = strokes_analyse(indices_known_creation);
    score_vals = [strokes_analyse(:).score];
   
    if ~isempty(strokes_analyse)
%         mask = ([strokes_analyse(:).assigned] >= last_stroke) & ( score_vals > 0.5 );
        mask = ([strokes_analyse(:).assigned] >= last_stroke) ;
        percentage_assigned_last = sum(mask)/length(strokes_analyse)*100;
    else
        percentage_assigned_last = NaN;
    end
    
end