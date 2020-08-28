function save_as_svg_sketch_visual_bad_strokes(strokes_topology, folder_save, filename, cam_param, set_removed, R)
%% Normilize pressure:
if ~isfield(strokes_topology(1), 'mean_pressure')
  for  i = 1:length(strokes_topology)
        strokes_topology(i).mean_pressure = mean(cat(1,strokes_topology(i).points2D(:).p));
  end
end


% pressure_max = max(cat(1, strokes_topology(:).mean_pressure));
% max_val = 0.75;
% if pressure_max < max_val 
%     scale_f = max_val / pressure_max;
%     for  i = 1:length(strokes_topology)
%         strokes_topology(i).mean_pressure = strokes_topology(i).mean_pressure*scale_f;
%     end
% end
%     

global sketch_height;
% global folder_save;
% 
% folder_save = fullfile(folder_save, 'lines_separation');
% if ~exist(folder_save, 'dir')
%    mkdir(folder_save); 
% end

filepath = fullfile(folder_save, filename);

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


for i = 1:length(strokes_topology)
   if ~strokes_topology(i).depth_assigned
       continue;
   end
   


   avg_pressure = strokes_topology(i).mean_pressure
   
   if avg_pressure < 1e-19
       fprintf('Zero pressure stroke %d\n', i);
       continue;       
   end
   
    points3D = strokes_topology(i).points3D(1,:)';
    
    point2D = cam_param.P*[R*points3D; 1.0];
    point2D = point2D(1:2)./point2D(3);
   
    str = sprintf('<path d="M %.5f %.5f ',...
                   point2D(1), ...
                   point2D(2));   
    fwrite(fid, str);
   
%    str = sprintf('L %.5f %.5f ',strokes_topology(i).primitive_geom(2),strokes_topology(i).primitive_geom(4));   
%    fwrite(fid, str);

    
    
   for j = 2:size(strokes_topology(i).points3D,1)
       points3D = strokes_topology(i).points3D(j,:)';
       point2D = cam_param.P*[R*points3D; 1.0];
       point2D = point2D(1:2)./point2D(3);
      str = sprintf('L %.5f %.5f ', point2D(1), ...
                   point2D(2));   
      fwrite(fid, str);
   end
   
   str = sprintf('"');   
   fwrite(fid, str);
   
   
   color = [0 0 0];
   
   if ismember(i, set_removed)
       color = [255 0 0];
   end
%   if ~IPad
    str = sprintf(' fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" id = "%d" />',...
                color(1),...
                color(2),...
                color(3),...
                pen_width*avg_pressure,...
                avg_pressure, ...
                i);   
%    else
%         str = sprintf(' fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" />',...
%                 color(1),...
%                 color(2),...
%                 color(3),...
%                 min(pen_width*avg_pressure, pen_width) + 1.0 ,...
%                 avg_pressure*3.3);  
%    end
   fwrite(fid, str);
end

str = sprintf('</svg>');
fwrite(fid, str);

fclose(fid);
end