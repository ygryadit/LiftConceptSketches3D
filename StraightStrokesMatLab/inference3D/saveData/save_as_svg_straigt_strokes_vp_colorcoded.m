function save_as_svg_straigt_strokes_vp_colorcoded(strokes_topology)

%% Normilise pressure:
pressure_max = max(cat(1, strokes_topology(:).mean_pressure));
max_val = 0.75;
if pressure_max < max_val 
    scale_f = max_val /pressure_max;
    for  i = 1:length(strokes_topology)
        strokes_topology(i).mean_pressure = strokes_topology(i).mean_pressure*scale_f;
    end
end
%%

global sketch_height;
global folder_save;
global last_added_stroke;
folder_save_svg = fullfile(folder_save, 'svg_files');
if ~exist(folder_save_svg, 'dir')
   mkdir(folder_save_svg); 
end



filepath = fullfile(folder_save_svg, sprintf('lines_vp.svg'));

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
% colors = uint8(colormap(parula(length(strokes)))*255);
colors = uint8(zeros(length(strokes_topology),3)*255);

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
       strokes_topology(i).points2D(1).x, ...
       strokes_topology(i).points2D(1).y);   
   fwrite(fid, str);
   
   for j = 2:length(strokes_topology(i).points2D)
      str = sprintf('L %.5f %.5f ',strokes_topology(i).points2D(j).x, strokes_topology(i).points2D(j).y);   
      fwrite(fid, str);
   end
   
   str = sprintf('"');   
   fwrite(fid, str);
   
   if strokes_topology(i).line_group == 1 
        color = uint8([1.0 0 0.5]*255);
   elseif strokes_topology(i).line_group == 2
        color = uint8([0 0.8 0.4]*255);   
   elseif strokes_topology(i).line_group == 3 
        color = uint8([0.4 0.3 1.0]*255);
   else
        color = uint8([0.8 0.8 0.0]*255);
   end
   
   str = sprintf(' fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" />',...
                color(1),...
                color(2),...
                color(3),...
                pen_width*avg_pressure,...
                avg_pressure);   
   fwrite(fid, str);
end

str = sprintf('</svg>');
fwrite(fid, str);

fclose(fid);
end