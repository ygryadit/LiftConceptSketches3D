function  percentage_low_score = computePecentageLowScoreAssignedLast(set_removed, strokes_topology_rough)
    
    strokes_analyse = strokes_topology_rough(set_removed);
   
    last_stroke = max([strokes_topology_rough(:).created]);
    
    
    indices_known_creation = [];
    for i = 1:length(strokes_analyse)
        if ~isempty(strokes_analyse(i).created)
            indices_known_creation(end+1) = i;
        end
    end
    
    strokes_analyse = strokes_analyse(indices_known_creation);
    
    confidence_vals = [strokes_analyse(:).confidence];
    score_vals = [strokes_analyse(:).score];
   
    if ~isempty(strokes_analyse)
        mask = ([strokes_analyse(:).assigned] >= last_stroke) & ( score_vals <=0.5 );
        percentage_low_score = sum(mask)/length(score_vals)*100;
    else
        percentage_low_score = NaN;
    end
    
end