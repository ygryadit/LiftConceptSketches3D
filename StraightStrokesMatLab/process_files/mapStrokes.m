function map = mapStrokes(stroke_topology, strokes)
map = NaN*ones(length(strokes),1);

ind2 = 1;

num_stroke_diff = abs(length(strokes) - length(stroke_topology)) + 1;
dist_thr = 10;

for i = 1:length(strokes)
    dist = Inf;
    
    ind2 = i;
    
    while (dist > dist_thr) & ((ind2 - i) < num_stroke_diff)
        if length(strokes(i).points) == length(stroke_topology(ind2).points2D)
                    
            dist = getDistPerPoint(strokes(i).points, stroke_topology(ind2).points2D);
        else
            %Check subsets
             dist = getDistPerPointSubset(strokes(i).points, stroke_topology(ind2).points2D);
        end

        if dist > dist_thr             
            ind2 = ind2 + 1;
        end
        
    end   
    if dist < dist_thr   
        map(i) = ind2;
    end
%     fprintf('i = %d, map = %d\n', i, map(i));
end

end

function dist = getDistPerPoint(points1, points2)
   dist1 = mean(sqrt(([points1(:).x] - [points2(:).x]).^2 + ([points1(:).y] - [points2(:).y]).^2));
   dist2 = mean(sqrt(([points1(end:-1:1).x] - [points2(:).x]).^2 + ...
                     ([points1(end:-1:1).y] - [points2(:).y]).^2));
   
   dist = min([dist1,dist2]);
end


function dist = getDistPerPointSubset(points1, points2)
    if length(points2) < length(points1)
        temp = points1;
        points1 = points2;
        points2 = temp;
    end
   dist = Inf;
   for j = 1:(length(points2) - length(points1) + 1)  
       dist_ = getDistPerPoint(points1, points2(j:(j+length(points1)-1)));
       dist = min([dist,dist_]);
   end

end