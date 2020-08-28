function save_as_svg_straigt_strokes_intersections_likely_v2(strokes, intersections, inds_intersections_corners, folder_save, sketch_height, name)


inds_others = setdiff(find(intersections.likely),inds_intersections_corners);    


pressure_max = max(cat(1, strokes(:).mean_pressure));
thr_max_pressure = 0.75;
if pressure_max < thr_max_pressure 
    scale_f = thr_max_pressure/pressure_max;
    for  i = 1:length(strokes)
        strokes(i).mean_pressure = strokes(i).mean_pressure*scale_f;
    end
end
    

filepath = fullfile(folder_save, sprintf(name));

pen_width = 3.0;
fid = fopen(filepath, 'w');

str = sprintf('<?xml version="1.0" encoding="utf-8" ?>\n'); 
fwrite(fid, str);
str = sprintf('<svg baseProfile="full" height="%d" version="1.1" viewBox="0,0,%d,%d" width="%d" xmlns="http://www.w3.org/2000/svg" xmlns:ev="http://www.w3.org/2001/xml-events" xmlns:xlink="http://www.w3.org/1999/xlink"><defs><style type="text/css"><![CDATA[\n', ...
    sketch_height,sketch_height,sketch_height,sketch_height);
fwrite(fid, str);
str = sprintf('\t.background { fill: white; }\n');
fwrite(fid, str);
str = sprintf('\t.line { stroke: firebrick; stroke-width: .1mm; }\n');
fwrite(fid, str);
str = sprintf('\t.blacksquare { fill: indigo; }\n');
fwrite(fid, str);
str = sprintf('\t.whitesquare { fill: white; }\n');
fwrite(fid, str);
str = sprintf(']]></style></defs>\n');
fwrite(fid, str);

  str = sprintf('<defs>\n');
  fwrite(fid, str);
  str = sprintf('<linearGradient id=''grad''>\n');
  fwrite(fid, str);
  str = sprintf('<stop stop-color=''red'' offset="0"/>\n');
  fwrite(fid, str);
  str = sprintf('<stop stop-color=''black'' offset="1"/>\n');
  fwrite(fid, str);
  str = sprintf('</linearGradient>\n');
  fwrite(fid, str);
  str = sprintf('</defs>\n');
  fwrite(fid, str);

%% Sketch:
for i = 1:length(strokes)   
   
   avg_pressure = strokes(i).mean_pressure;
   
   if avg_pressure < 1e-19
       fprintf('Zero pressure stroke %d\n', i);
       continue;       
   end
   
   str = sprintf('<path d="M %.5f %.5f ',...
       strokes(i).points2D(1).x, ...
       strokes(i).points2D(1).y);   
   fwrite(fid, str);
   
   for j = 2:length(strokes(i).points2D)
      str = sprintf('L %.5f %.5f ',strokes(i).points2D(j).x, strokes(i).points2D(j).y);   
      fwrite(fid, str);
   end
   
   str = sprintf('"');   
   fwrite(fid, str);

   color = [0 0 0];
   str = sprintf(' fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" />',...
                color(1),...
                color(2),...
                color(3),...
                pen_width*avg_pressure,...
                avg_pressure);   
   fwrite(fid, str);
end


%% Intersections:


% circle_r = 

for i = 1:length(intersections.likely)
%     if ~intersections(i).likely
%         continue;
%     end
        
    if ismember(i, inds_intersections_corners)
        color = [255, 153, 255];
    elseif ismember(i, inds_others)
        color = [33, 150, 243]; 
    else
        continue;
    end
    
    if sum([strokes(intersections.strokes_indices(i,:)).primitive_type] == 0) ~=2
        continue;
    end
    
    str = sprintf(' <circle cx="%.3f" cy="%.3f" r="%.3f" stroke="black" fill-opacity="0.5" fill="rgb(%d, %d, %d)" />',...
                intersections.coordinates2D(i,1),...
                intersections.coordinates2D(i,2),...
                intersections.accuracy_radius(i),...              
                color(1),...
                color(2),...
                color(3));   
   fwrite(fid, str);           
end




str = sprintf('</svg>');
fwrite(fid, str);

fclose(fid);
end

function coord = findClosestTwoPointsStroke(strokePoints, point)
    strokePoints = [cat(1,strokePoints.x) cat(1,strokePoints.y)];
    num_points = size(strokePoints,1);
    point = repmat(point, num_points, 1);
    distances = sqrt(sum((strokePoints - point).^2,2));
    [~, inds] = sort(distances);
    coord(1,:) = strokePoints(inds(1),:);
    point = point(1,:);
    dir1 = coord(1,:) - point;
    step = min(20, norm(dir1));
    coord(1,:) = point + step*dir1./norm(dir1);
    
    coord(2,:) = strokePoints(inds(2),:);
    dir2 = coord(2,:) - point;
    step = min(20, norm(dir2));
    coord(2,:) = point + step*dir2./norm(dir2);
end

function color = agrement2color(value)
colors = parula(256);
% color = uint8(([1.0 1.0 1.0] - colors(uint8((1-value)*255)+1,:))*256);
% color = uint8((colors(uint8(value*255)+1,:))*256);
color = uint8((colors(uint8((1.0 - value)*255)+1,:))*256);
% color = [255 0 0];
end