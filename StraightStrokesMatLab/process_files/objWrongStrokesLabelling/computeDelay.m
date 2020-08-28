function delays = computeDelay(set_removed, strokes_topology_rough)
    
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
        delays = [strokes_analyse(:).assigned] - [strokes_analyse(:).created];       
    else
        delays = NaN;
    end
    delays = delays';
end