function save_as_svg_full_sketch(strokes_topology, intersections, folder_save, sketch_height)
global IPad;
%% Normilize pressure:
if ~isfield(strokes_topology(1), 'mean_pressure')
  for  i = 1:length(strokes_topology)
        strokes_topology(i).mean_pressure = mean(cat(1,strokes_topology(i).points2D(:).p));
  end
end


pressure_max = max(cat(1, strokes_topology(:).mean_pressure));
max_val = 0.75;
if pressure_max < max_val 
    scale_f = max_val /pressure_max;
    for  i = 1:length(strokes_topology)
        strokes_topology(i).mean_pressure = strokes_topology(i).mean_pressure*scale_f;
    end
end
    


% global folder_save;
% 
% folder_save = fullfile(folder_save, 'lines_separation');
% if ~exist(folder_save, 'dir')
%    mkdir(folder_save); 
% end

filepath = fullfile(folder_save, 'intersections_active.svg');

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
colors = uint8(colormap(parula(length(strokes_topology)))*255);


for i = 1:length(strokes_topology)
%    if ~ismember(strokesType(i), [contour, key])         
%        continue;
%    end
   


   avg_pressure = strokes_topology(i).mean_pressure;
   
   if avg_pressure < 1e-19
       fprintf('Zero pressure stroke %d\n', i);
       continue;       
   end
   
   str = sprintf('<path d="M %.5f %.5f ',...
                   strokes_topology(i).points2D(1).x, ...
                   strokes_topology(i).points2D(1).y);   
   fwrite(fid, str);
   
%    str = sprintf('L %.5f %.5f ',strokes_topology(i).primitive_geom(2),strokes_topology(i).primitive_geom(4));   
%    fwrite(fid, str);
   for j = 2:length(strokes_topology(i).points2D)
      str = sprintf('L %.5f %.5f ',strokes_topology(i).points2D(j).x, strokes_topology(i).points2D(j).y);   
      fwrite(fid, str);
   end
   
   str = sprintf('"');   
   fwrite(fid, str);
   
   
   color = [0 0 0];
   
%    color = colors(i,:);
   

     str = sprintf(' fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" id = "%d" /> \n',...
                color(1),...
                color(2),...
                color(3),...
                pen_width*avg_pressure,...
                avg_pressure, ...
                i);   

   fwrite(fid, str);
end

%% Plot active intersections
radius_circle = 3;
for i = 1:length(intersections)
   if isempty(intersections(i).is_active)
       continue;
   end
   
   if intersections(i).is_active
       color = [0 255 0];
   else
       color = [255 0  0];
   end
      
    str = sprintf(' <circle cx="%.3f" cy="%.3f" r="%.3f" stroke="black" fill-opacity="0.5" fill="rgb(%d, %d, %d)" />',...
                intersections(i).coordinates2D(1),...
                intersections(i).coordinates2D(2),...
                radius_circle,...                
                color(1),...
                color(2),...
                color(3));   
   fwrite(fid, str);  
end


str = sprintf('</svg>');
fwrite(fid, str);

fclose(fid);



end

function plot_radius(fid, strokes_topology, i, color)
str = sprintf('<path d="M %.5f %.5f ',...
                   strokes_topology(i).points2D(1).x, ...
                   strokes_topology(i).points2D(1).y);   
fwrite(fid, str);
for j = 2:length(strokes_topology(i).points2D)
      str = sprintf('L %.5f %.5f ',strokes_topology(i).points2D(j).x, strokes_topology(i).points2D(j).y);   
      fwrite(fid, str);
end


str = sprintf('" fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" stroke-linecap="round" id = "%d" /> \n',...
            color(1),...
            color(2),...
            color(3),...
            5,...
            0.5, ...
            i);   

fwrite(fid, str);

end