function [strokes_topology, intersections, cam_param] = ...
            scaleAndCenterTheObject(strokes_topology, intersections, cam_param, img)
% folder_results = 'D:\Projects\WiresProject\Results\output_v11';
% if ~exist('view', 'var')
%     view = 'view1';
% end
% filepath = fullfile(folder_results,prof,obj,view,'Professional2_vacuum_cleaner_full.json');
% [strokes_topology, intersections] = readReconstructionJson(filepath);
% 
% global ZUP
% if ~ZUP
%     signz = -1;
%     ZUP = true;
% else
%     signz = 1; 
% end

%% Center:
[strokes_topology, intersections, cam_param] = ...
                centerObject(strokes_topology, intersections, cam_param, img);


[strokes_topology, intersections, cam_param] = ...
    scaleObject(strokes_topology, intersections, cam_param, img);

% strokes_topology = clipOutsideBB(strokes_topology, intersections, cam_param);


pressure_max = max(cat(1, strokes_topology(:).mean_pressure));

thr_max_pressure = 0.75;
if pressure_max < thr_max_pressure 
    scale_f = thr_max_pressure/pressure_max;
    for  i = 1:length(strokes_topology)
        strokes_topology(i).mean_pressure = strokes_topology(i).mean_pressure*scale_f;
    end
end

end



