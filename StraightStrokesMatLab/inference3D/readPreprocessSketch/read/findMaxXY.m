function dim = findMaxXY(sketch)

dim_max  = 0;
dim_min = Inf;

for i = 1:length(sketch.strokes)
    vals = [sketch.strokes(i).points(:).x sketch.strokes(i).points(:).y];
    dim_max = max([ dim_max vals]);
    dim_min = min([ dim_min vals]);
    
end

dim = dim_max+dim_min;


end