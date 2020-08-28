function map = toMerged(strokes_topology)
map = NaN*ones(length(strokes_topology),1);

for i = 1:length(strokes_topology)
    if isempty(strokes_topology(i).merged_with)
        map(i) = i;
    else
        map(i) = min(strokes_topology(i).merged_with);
    end
end

end