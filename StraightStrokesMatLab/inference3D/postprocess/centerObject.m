function [strokes_topology, intersections, cam_param] = ...
                centerObject(strokes_topology, intersections, cam_param, img)

img = ones(512,512,3);
reproject3Dto2D(img, cam_param, strokes_topology,intersections,1,'k');

bound_box = bbOFIntersections(intersections);
center = mean(bound_box,2)';




%% Stroke_topology:
for i = 1:length(strokes_topology)
    if ~isempty(strokes_topology(i).points3D)
        strokes_topology(i).points3D  = ...
            strokes_topology(i).points3D - ...
            center;
    end
    
    if isfield(strokes_topology(i), 'primitive_geom_3D') & ~isempty(strokes_topology(i).primitive_geom_3D)
        strokes_topology(i).primitive_geom_3D = ...
            strokes_topology(i).primitive_geom_3D - ...
            center;
    end
end

%% Intersections:


for i = 1:length(intersections)
    if isnan(intersections(i).coordinates3D)
        continue;
    end
    
    intersections(i).coordinates3D = intersections(i).coordinates3D - center;    
    
%     p2D = cam_param.P*[intersections(i).coordinates3D'; 1.0];
% 
%     p2D = p2D(1:2)./p2D(3);
%             
%     plot(p2D(1),p2D(2), '*');
end

% reproject3Dto2D(img, cam_param, strokes_topology,intersections,1,'r');
%% Camera parameters

cam_param.K
cam_param.R
cam_param.C

%recalculate rotation matrix:
up = [0,0,1];
% cam_pos = cam_param.C;
% focal_point = cam_param.C-cam_param.view_dir';
% cam_param.R = rotationMatrixFromView(cam_pos', focal_point', up);


% reproject3Dto2D(img, cam_param, strokes_topology,intersections,2);
cam_param.C = cam_param.C - center';
cam_pos = cam_param.C;
cam_pos = reshape(cam_pos, 1, 3);
view_dir = reshape(cam_param.view_dir, 1, 3);
focal_point = cam_pos+view_dir;
cam_param.R
% cam_param.R = rotationMatrixFromView(cam_pos, focal_point, up);
cam_param.R = rotationMatrixFromView(cam_pos, focal_point, cam_param.R(2,:));
cam_param.R
cam_param.P =  cam_param.K *[ cam_param.R -cam_param.R*cam_param.C];
cam_param.t = -cam_param.R*cam_param.C;


reproject3Dto2D(img, cam_param, strokes_topology,intersections,1, 'r:');

end