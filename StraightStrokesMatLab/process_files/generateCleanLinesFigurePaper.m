%% read felix's data:
folder = "C:\Users\yulia\Research\DesignSketch3D\paper_last\images_draft";
filename = "designer2_printer_02_bestScore_full_with_clustering.json";

[ strokes_topology_curves, intersections, cam_param] = ...
            readReconstructionJsonFelix( fullfile(folder, filename));   
        
%% read the original file with clean curves:    
folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\_400';
designer = 'designer2';
object = 'printer_02';
view = 'view1';
filename = 'designer2_printer_02_bestScore_full.json';

[ strokes_topology_rough, ~, ~] = ...
            readReconstructionJson( fullfile(folder_files, designer, object, view, filename)); 
        
%% Save clean image:
% read canvas height:
folder_original_file = 'C:\Users\yulia\Research\Data\IndustrialDesignSketches\sketches_json_first_viewpoint';
filepath_sketch_json = fullfile(folder_original_file, designer, object, 'view1_concept.json');
sketch = readSketchJson(filepath_sketch_json);
sketch_height = sketch.canvas.width;
clear('sketch');

filepath_save = fullfile(folder, 'designer2_printer_02_bestScore_full_input.svg');

%% Header
pen_width = 1.5;
fid = fopen(filepath_save, 'w');

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



for i = 1:length(strokes_topology_rough)
%    if ~strokes_topology_rough(i).depth_assigned
%        continue;
%    end
   
   % consider only straigth strokes:
%    if strokes_topology_rough(i).primitive_type ~= 0 && strokes_topology_rough(i).primitive_type ~= -2
%        continue;
%    end
%    
   avg_pressure = strokes_topology_rough(i).mean_pressure;
 
   if avg_pressure < 1e-19
       fprintf('Zero pressure stroke %d\n', i);
       continue;       
   end
   
%    if  ~isempty(strokes_topology_rough(i).merged_with) && (min(strokes_topology_rough(i).merged_with) ~= i)
%       continue;      
%    end
   
   if length(strokes_topology_rough(i).points2D) < 2
        continue;
   end
   point2D(1) = strokes_topology_rough(i).points2D(1).x;
   point2D(2) = strokes_topology_rough(i).points2D(1).y;
   
   str = sprintf('<path d="M %.5f %.5f ',...
                   point2D(1), ...
                   point2D(2));   
   fwrite(fid, str);
   
   
   for j = 2:size(strokes_topology_rough(i).points3D,1)
      point2D(1) = strokes_topology_rough(i).points2D(j).x;
      point2D(2) = strokes_topology_rough(i).points2D(j).y;
   
      str = sprintf('L %.5f %.5f ', point2D(1), ...
                   point2D(2));   
      fwrite(fid, str);
   end
   
   str = sprintf('"');   
   fwrite(fid, str);
   
   
   color = [0 0 0];

   str = sprintf(' fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" id = "%d" />',...
                color(1),...
                color(2),...
                color(3),...
                pen_width*avg_pressure,...
                avg_pressure, ...
                i);   

   fwrite(fid, str);
end

% for i = 1:length(strokes_topology_curves)
% %    if ~strokes_topology_curves(i).depth_assigned
% %        continue;
% %    end
%    
%    % consider only straigth strokes:
%    if strokes_topology_curves(i).primitive_type ~= 1
%        continue;
%    end
%    
%    avg_pressure = strokes_topology_curves(i).mean_pressure;
% 
%   
%     points3D = strokes_topology_curves(i).points3D(1,:)';
%     point2D = cam_param.P*[points3D; 1.0];
%     point2D = point2D(1:2)./point2D(3);
%    
%    str = sprintf('<path d="M %.5f %.5f ',...
%                    point2D(1), ...
%                    point2D(2));   
%    fwrite(fid, str);
%    
% 
%    for j = 2:size(strokes_topology_curves(i).points3D,1)
%        points3D = strokes_topology_curves(i).points3D(j,:)';
%        point2D = cam_param.P*[points3D; 1.0];
%        point2D = point2D(1:2)./point2D(3);
%    
%       str = sprintf('L %.5f %.5f ', point2D(1), ...
%                    point2D(2));   
%       fwrite(fid, str);
%    end
%    
%    str = sprintf('"');   
%    fwrite(fid, str);
% 
%    color = [0 0 0];
% 
%    str = sprintf(' fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" id = "%d" />',...
%                 color(1),...
%                 color(2),...
%                 color(3),...
%                 pen_width*avg_pressure,...
%                 avg_pressure, ...
%                 i);   
% 
%    fwrite(fid, str);
% end


str = sprintf('</svg>');
fwrite(fid, str);
fclose(fid);


