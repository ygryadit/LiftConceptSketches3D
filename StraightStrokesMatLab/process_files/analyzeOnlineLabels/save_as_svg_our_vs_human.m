function save_as_svg_our_vs_human(strokes,intersections_coordinates,solution_labeling, solution_our, folder_save, sketch_height, name)





% [strokes] = computeAverageStrokesSpeed(strokes);
% [accuracy_radiuses] = computeMergeThresholdStrokes(cat(1,strokes(:).speed));
radius_circle = 2.5;
 
num_strokes = round(length(strokes)); 
   

pressure_max = max(cat(1, strokes(:).mean_pressure));
thr_max_pressure = 0.6;
if pressure_max < thr_max_pressure 
    scale_f = thr_max_pressure/pressure_max;
    for  i = 1:length(strokes)
        strokes(i).mean_pressure = strokes(i).mean_pressure*scale_f;
    end
end
    

filepath = fullfile(folder_save, sprintf(name));

pen_width = 1.5;
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
%    
%    if strokes(i).line_group == 1 
%         color = uint8([1.0 0 0.5]*255);
%    elseif strokes(i).line_group == 2
%         color = uint8([0 0.8 0.4]*255);   
%    elseif strokes(i).line_group == 3 
%         color = uint8([0.4 0.3 1.0]*255);
%    else
%         color = uint8([0.8 0.8 0.0]*255);
%    end
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



agrement = (solution_labeling == solution_our');

% circle_r = 

for i = 1:size(intersections_coordinates,1)
    color = agrement2color(agrement(i));
    str = sprintf(' <circle cx="%.3f" cy="%.3f" r="%.3f" stroke="black" fill-opacity="0.5" fill="rgb(%d, %d, %d)" />',...
                intersections_coordinates(i,1),...
                intersections_coordinates(i,2),...
                radius_circle,...                
                color(1),...
                color(2),...
                color(3));   
   fwrite(fid, str);         
%    strk_indices = intersections.strokes_indices(i,:);
%    
%    
%    coord = findClosestTwoPointsStroke( strokes(strk_indices(1)).points2D,...
%                                     intersections_coordinates(i,:)); 
%    
% %   str = sprintf('<defs>\n');
% %   fwrite(fid, str);
% %   str = sprintf('<linearGradient id=''grad%d''  x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f">\n', i, ...
% %                      intersections_coordinates(i,1),...
% %                      intersections_coordinates(i,2),...
% %                      coord(1,1),...
% %                      coord(1,2));
% %   fwrite(fid, str);
% %   str = sprintf('<stop stop-color="rgb(%d, %d, %d)" offset="0"/>\n',...
% %                 color(1),...
% %                 color(2),...
% %                 color(3));
% %   fwrite(fid, str);
% %   str = sprintf('<stop stop-color=''black'' offset="1"/>\n');
% %   fwrite(fid, str);
% %   str = sprintf('</linearGradient>\n');
% %   fwrite(fid, str);
% %   str = sprintf('</defs>\n');
% %   fwrite(fid, str);
% %   
% 
%     p = 0.9;
%     sc = 2.6;
%    avg_pressure = strokes(strk_indices(1)).mean_pressure;
%    str = sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" />',...
%                 intersections_coordinates(i,1),...
%                 intersections_coordinates(i,2),...
%                 coord(1,1),...
%                 coord(1,2),...
%                 color(1),...
%                 color(2),...
%                 color(3),...
%                 pen_width*avg_pressure*sc,...
%                 avg_pressure);
%    fwrite(fid, str);
%    
%    
%    
%    str = sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" />',...
%                 intersections_coordinates(i,1),...
%                 intersections_coordinates(i,2),...
%                 coord(2,1),...
%                 coord(2,2),...
%                 color(1),...
%                 color(2),...
%                 color(3),...
%                 pen_width*avg_pressure*sc,...
%                 avg_pressure);
%    fwrite(fid, str);
%    
%    
%    coord = findClosestTwoPointsStroke( strokes(strk_indices(2)).points2D,...
%                                     intersections_coordinates(i,:));
%    avg_pressure = strokes(strk_indices(2)).mean_pressure;
%    str = sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" />',...
%                 intersections_coordinates(i,1),...
%                 intersections_coordinates(i,2),...
%                 coord(1,1),...
%                 coord(1,2),...
%                 color(1),...
%                 color(2),...
%                 color(3),...
%                 pen_width*avg_pressure*sc,...
%                 avg_pressure);
%    fwrite(fid, str);
%    
%    
%    
%    str = sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" />',...
%                 intersections_coordinates(i,1),...
%                 intersections_coordinates(i,2),...
%                 coord(2,1),...
%                 coord(2,2),...
%                 color(1),...
%                 color(2),...
%                 color(3),...
%                 pen_width*avg_pressure*sc,...
%                 avg_pressure);
%    fwrite(fid, str);
   
   
   
   
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

% color = uint8(([1.0 1.0 1.0] - colors(uint8((1-value)*255)+1,:))*256);
% color = uint8((colors(uint8(value*255)+1,:))*256);
if value
    color =[0 255 0];
else
    color =[255 0 0];
end
% color = uint8((colors(uint8((1.0 - value)*255)+1,:))*256);
% color = [255 0 0];
end