function save_as_svg_primitives(strokes_topology, filename)

%% Normilize pressure:
pressure_max = max(cat(1, strokes_topology(:).mean_pressure));
max_val = 0.9;
if pressure_max < max_val 
    scale_f = max_val /pressure_max;
    for  i = 1:length(strokes_topology)
        strokes_topology(i).mean_pressure = strokes_topology(i).mean_pressure*scale_f;
    end
end
    

global sketch_height;
global folder_save;

folder_save_svg = fullfile(folder_save, 'merged_svg');



if ~exist(folder_save_svg, 'dir')
   mkdir(folder_save_svg); 
end

filepath = fullfile(folder_save_svg, filename);

global fid;
fprintf(fid, 'Saved to %s\n', filepath);

pen_width = 1.5;
fids = fopen(filepath, 'w');

str = sprintf('<?xml version="1.0" encoding="utf-8" ?>\n'); 
fwrite(fids, str);
str = sprintf('<svg baseProfile="full" height="%d" version="1.1" viewBox="0,0,%d,%d" width="%d" xmlns="http://www.w3.org/2000/svg" xmlns:ev="http://www.w3.org/2001/xml-events" xmlns:xlink="http://www.w3.org/1999/xlink"><defs><style type="text/css"><![CDATA[\n', ...
    sketch_height,sketch_height,sketch_height,sketch_height);
fwrite(fids, str);
str = sprintf('\t.background { fill: white; }\n');
fwrite(fids, str);
str = sprintf('\t.line { stroke: firebrick; stroke-width: .1mm; }\n');
fwrite(fids, str);
str = sprintf('\t.blacksquare { fill: indigo; }\n');
fwrite(fids, str);
str = sprintf('\t.whitesquare { fill: white; }\n');
fwrite(fids, str);
str = sprintf(']]></style></defs>\n');
fwrite(fids, str);
% colors = uint8(colormap(parula(length(strokes)))*255);


for i = 1:length(strokes_topology)
%    if ~ismember(strokesType(i), [contour, key])         
%        continue;
%    end
   
   if strokes_topology(i).primitive_type ~= 0
       continue;
   end

   avg_pressure = strokes_topology(i).mean_pressure;
   
   if avg_pressure < 1e-19
       fprintf('Zero pressure stroke %d\n', i);
       continue;       
   end
   
   str = sprintf('<path d="M %.5f %.5f ',...
                   strokes_topology(i).primitive_geom(1), ...
                   strokes_topology(i).primitive_geom(3));   
   fwrite(fids, str);
   
   str = sprintf('L %.5f %.5f ',strokes_topology(i).primitive_geom(2),strokes_topology(i).primitive_geom(4));   
   fwrite(fids, str);
%    for j = 2:length(strokes_topology(i).points2D)
%       str = sprintf('L %.5f %.5f ',strokes_topology(i).points2D(j).x, strokes_topology(i).points2D(j).y);   
%       fwrite(fids, str);
%    end
   
   str = sprintf('"');   
   fwrite(fids, str);
   
   
   color = [0 0 0];
   
   
   str = sprintf(' fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" />',...
                color(1),...
                color(2),...
                color(3),...
                pen_width*avg_pressure,...
                avg_pressure);   
   fwrite(fids, str);
end

str = sprintf('</svg>');
fwrite(fids, str);

fclose(fids);
end