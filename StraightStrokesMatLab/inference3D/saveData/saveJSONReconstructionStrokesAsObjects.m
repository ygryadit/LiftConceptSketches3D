function saveJSONReconstruction(strokes_topology, ...
                                intersections,...
                                cam_param,...
                                folder_save, ...
                                fname)
global designer;
global object_name;
global last_added_stroke;

if ~exist('fname', 'var')
    fname = 'confident';
    
end

global SAVE_SVGS;
if ~SAVE_SVGS
    return;
end

global ZUP
if ~ZUP
    signz = -1;
else
    signz = 1; 
end
fprintf('signz = %d\n',signz);

% for i = 1:length(strokes_topology)
% %     strokes_topology(i).points3D(:,3) = -strokes_topology(i).points3D(:,3)*signz;
%     strokes_topology(i).points3D = strokes_topology(i).points3D*;
%     strokes_topology(i).primitive_geom_3D = strokes_topology(i).primitive_geom_3D;
% end
% 
% for i = 1:length(intersections)
% %     intersections(i).coordinates3D(:,3) = -intersections(i).coordinates3D(:,3)*signz;
%     intersections(i).coordinates3D = intersections(i).coordinates3D*signz;
% end

object.designer     = designer;
object.object_name  = object_name;
object.cam_param    = cam_param;
object.strokes_topology = strokes_topology;
object.intersections = intersections;
object.signz = signz;
object_str = jsonencode(object);

fullfile(folder_save, ...
                     sprintf('%s_%s_%s_full.json', ...
                     designer,...
                     object_name,...
                     fname))
                 
fid_ = fopen(fullfile(folder_save, ...
                     sprintf('%s_%s_%s_full.json', ...
                     designer,...
                     object_name,...
                     fname)), ...
                     'w');
fwrite(fid_, object_str);
fclose(fid_); 

% fid = fopen(fullfile(folder_save, sprintf('%s_%s_%d_full.json',designer, object_name, last_added_stroke)), 'w');
% fwrite(fid, object_str);
% fclose(fid); 

clear('object');
object.designer  = designer;
object.object_name  = object_name;

% %% Light
% 
% 
% for i = 1:length(strokes_topology)
% strokes(i).points2D        = strokes_topology(i).points2D;
% strokes(i).points3D        = strokes_topology(i).points3D;
% strokes(i).mean_pressure   = strokes_topology(i).mean_pressure;
% if strokes_topology(i).primitive_type==0
%     strokes(i).line_group      = strokes_topology(i).line_group;
% else
%     strokes(i).line_group      = -1;
% end
% strokes(i).depth_assigned  = strokes_topology(i).depth_assigned;
% end
% 
% intersections_ = intersections;
% clear('intersections');
% for i = 1:length(intersections_)
%     intersections(i).coordinates2D  = intersections_(i).coordinates2D;
%     intersections(i).coordinates3D  = intersections_(i).coordinates3D;
%     intersections(i).collinear      = intersections_(i).collinear;
%     intersections(i).is_active      = intersections_(i).is_active;
%     intersections(i).strokes_indices  = intersections_(i).strokes_indices;
% end
% 
% object.intersections = intersections;
% object.strokes = strokes;
% object.cam_param = cam_param;
% object.cam_param = cam_param;
% 
% object_str = jsonencode(object);
% object_str = processJSONstr(object_str);
% 
% fid = fopen(fullfile(folder_save, [designer '_' object_name '_' fname '.json']), 'w');
% fwrite(fid, object_str);
% fclose(fid); 
% 
% 
% 
% %% Different order of coordinates 
% for i = 1:length(strokes)
%     for j = 1:size(strokes(i).points3D,1)
%         strokes(i).points3D(j,:)   =  strokes(i).points3D(j, [1,3,2]);
%         strokes(i).points3D(j,2)     = -strokes(i).points3D(j,2);
%     end
% end
% 
% for i = 1:length(intersections)   
%     intersections(i).coordinates3D(:)  = intersections(i).coordinates3D([1,3,2]);
%     intersections(i).coordinates3D(2) = -intersections(i).coordinates3D(2);
% end
% object.strokes = strokes;
% object.intersections = intersections;
% object_str           = jsonencode(object);
% object_str           = processJSONstr(object_str);
% 
% fid = fopen(fullfile(folder_save, [designer '_' object_name '_' fname '_y_up.json']), 'w');
% fwrite(fid, object_str);
% fclose(fid); 


%% Save SVG color coding:

filename = sprintf('%s_%s_%s_strokes_scores.svg', ...
                     designer,...
                     object_name,...
                     fname);
                 
save_as_svg_confidence_stroke(strokes_topology, folder_save, filename)
end

function str = processJSONstr(str)

str = strrep(str, ',', sprintf(',\r'));
str = strrep(str, '[{', sprintf('[\r{\r'));
str = strrep(str, '}]', sprintf('\r}\r]'));

end