function plotMapping(stroke_topology, strokes, map)
   
    figure(1);
for i = 1:length(strokes)

    hold off;
    if ~isnan(map(i))
        plot([strokes(i).points(:).x], [strokes(i).points(:).y],'r')    
        hold on;
        plot([stroke_topology(map(i)).points2D(:).x], [stroke_topology(map(i)).points2D(:).y],'b:')
        axis equal;
        axis ij;
    end
end

end