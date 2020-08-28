function [strokes_topology, intersections, cam_param] = ...
    scaleObject(strokes_topology, intersections, cam_param, img)

strokes_assigned = strokes_topology(cat(1, strokes_topology(:).depth_assigned));
plotStrokesTopology(3, strokes_assigned);
plot3(cam_param.C(1), cam_param.C(2), cam_param.C(3), '*r');
plot3([cam_param.C(1), cam_param.C(1)+cam_param.view_dir(1)], [cam_param.C(2), cam_param.C(2)+cam_param.view_dir(2)], [cam_param.C(3), cam_param.C(3)+cam_param.view_dir(3)], '-');

bound_box = bbOFIntersections(intersections);

dims = bound_box(:,2) - bound_box(:,1);
max_dim = max(dims);

% global ZUP
% if ~ZUP
%     signz = -1;
%     ZUP = true;
% else
%     signz = 1; 
% end

scale_factor = 0.4./max_dim;


%% Stroke_topology:
for i = 1:length(strokes_topology)
    if ~isempty(strokes_topology(i).points3D)
        strokes_topology(i).points3D  = ...
            strokes_topology(i).points3D*scale_factor;
    end
    
    if isfield(strokes_topology(i), 'primitive_geom_3D') & ~isempty(strokes_topology(i).primitive_geom_3D)
        strokes_topology(i).primitive_geom_3D = ...
            strokes_topology(i).primitive_geom_3D*scale_factor;
    end
end

%% Intersections:


for i = 1:length(intersections)
    if isnan(intersections(i).coordinates3D)
        continue;
    end
    
    intersections(i).coordinates3D = intersections(i).coordinates3D*scale_factor;    
    
end

%% Camera parameters
up = [0,0,1];

cam_param.C = cam_param.C*scale_factor;
cam_pos = cam_param.C;
cam_pos = reshape(cam_pos, 1, 3);
view_dir = reshape(cam_param.view_dir, 1, 3);
focal_point = cam_pos+view_dir;
cam_param.R
cam_param.R = rotationMatrixFromView(cam_pos, focal_point, cam_param.R(2,:));
cam_param.R
cam_param.P =  cam_param.K *[ cam_param.R -cam_param.R*cam_param.C];
cam_param.t = -cam_param.R*cam_param.C;


reproject3Dto2D(img, cam_param, strokes_topology,intersections,1, 'b.-');

strokes_assigned = strokes_topology(cat(1, strokes_topology(:).depth_assigned));
plotStrokesTopology(4, strokes_assigned);
plot3(cam_param.C(1), cam_param.C(2), cam_param.C(3), '*r');
plot3([cam_param.C(1), cam_param.C(1)+cam_param.view_dir(1)], [cam_param.C(2), cam_param.C(2)+cam_param.view_dir(2)], [cam_param.C(3), cam_param.C(3)+cam_param.view_dir(3)], '-');
end