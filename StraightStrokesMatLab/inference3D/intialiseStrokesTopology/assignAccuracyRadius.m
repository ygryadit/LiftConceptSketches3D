function strokes_topology = assignAccuracyRadius(strokes_topology, accuracy_radiuses)

    for i = 1:length(strokes_topology)        
        strokes_topology(i).accuracy_radius = accuracy_radiuses(i);
    end

end