function lines_groups = linesGroupsOrthographicCamera(strokes_topology)
 
    
   mask_lines          = cat(1, strokes_topology(:).primitive_type) == 0;
   lines_coordinates2D = cat(1,strokes_topology(mask_lines).primitive_geom);
   dirs = lines_coordinates2D(:,[2,4]) - lines_coordinates2D(:,[1,3]);
   
   %% Find sets of parallel lines:
   dirs = dirs';
   dirs = dirs./vecnorm(dirs);
   dirs = dirs';
   D2 = pdist(dirs,'cosine');
   Z = linkage(D2);
   T = cluster(Z,'cutoff',1.1);
   showCluster(T, lines_coordinates2D)
end


function dist = angilarDistance(d1, d2)
    dist = 1 - abs(dot(d1, d2));
end
 
function showCluster(T, lines_coordinates2D)
   global filepath_sketch_img;    
    img = readSketchImg(filepath_sketch_img);
colors = uint8(hsv(max(T))*255);
figure(1);
hold off;
imshow(img);
hold on;

for i = 1:max(T)
   inds = find(T == i);
   for j = inds'
      plot(lines_coordinates2D(j,[1,2]),lines_coordinates2D(j,[3,4]), ...
          'Color', colors(i,:));
   end
end

end