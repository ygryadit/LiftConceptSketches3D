function viewsForDesigner()
folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';

folder = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\curves17052020\results\results_medium';
global folder_save;
global sketch_height;
sketch_height = 512;

files = dir(folder);
files = files(~[files(:).isdir]);
files = {files.name};



foler_save_base = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\paper_viewpoints_medium\';

for i = 1:length(files)
   close all
   filename = files{i};
   filename_ = strrep(filename, '.json', '');
   members = split(filename_,'_');
   designer = members{1};
   
   object_parts = members(2:(end-2));
   object = object_parts{1};
   for j = 2:length(object_parts)
      object = sprintf('%s_%s', object, object_parts{j});
   end
   
   view = 'view1';
   
   [success, strokes_topology, intersections, cam_param] = ...
                combineRead(folder, filename,...
                            folder_files,designer, object, view);
                        
% viewsObject('designer2_guitar_01', foler_save_base, folder)
% viewsObject('Prof2task2_cabinet_01', foler_save_base, folder)
    close all
    if success
        viewsObject(filename_, foler_save_base, folder, strokes_topology, intersections, cam_param);
    end
end

end

function [success, strokes_topology, intersections, cam_param] = ...
                combineRead(folder, filename,...
                            folder_files,designer, object, view)
   success = true;
   try
       [ strokes_topology, intersections, cam_param] = ...
            readReconstructionJsonFelix( fullfile(folder, filename));   
   catch
       success = false;
       return;
   end
   
   try
   [ strokes_topology_, ~, ~] = ...
            readReconstructionJson( fullfile(folder_files,designer, object, view, filename));   
   catch
       success = false;
       return;
   end
        
    for ii = 1:length(strokes_topology)
        if strokes_topology_(ii).primitive_type ~= 0
            continue;
        end

        strokes_topology(ii).points3D = strokes_topology_(ii).points3D;
        if ~isempty(strokes_topology(ii).points3D)
            strokes_topology(ii).primitive_type = 0;
        end
    end
    
end



function viewsObject(filename, foler_save_base, folder, strokes_topology, intersections, cam_param)
    global folder_save;
   
    [strokes_topology, intersections, cam_param] = ...
        scaleAndCenterTheObject(strokes_topology, intersections, cam_param, []);
    folder_save = fullfile(foler_save_base, filename);

    
    cam_param.view_dir = cam_param.view_dir./norm(cam_param.view_dir);
    vec = cross(cam_param.view_dir, [0 0 1]');
    up = cross(cam_param.view_dir,vec);
    vec = cross(cam_param.view_dir, up);
    
    alpha = lsqnonlin(@distCenter, norm(cam_param.view_dir+cam_param.C));
    focal_point = cam_param.C + alpha * cam_param.view_dir;

    l = norm(cam_param.C - focal_point);
    start = cam_param.C - up*l*0.4;

    frame_num = 0;
    angle = 0;
    scaleX = 1.0;
    scaleY = 1.0;
    scaleZ = 1.0;
    
    
   
    
    
    
  
    for i = 0:45:315
        frame_num = frame_num+1;
        R=AxelRot(i,cam_param.view_dir, cam_param.C);
        s = R*[start; 1];
        figure(4);
        plot3(s(1), s(2), s(3), '*');
        cam_pos = s(1:3);
        folder_name = 'rotate';
        applyTransformSaveSVGs(strokes_topology, angle, scaleX, scaleY, scaleZ, focal_point, cam_pos, cam_param, frame_num, folder_name, intersections)
    end
    
%     frame_num = 0;
%     for i = 0:45:360
%         frame_num = frame_num+1;
%         R=AxelRot(i,cam_param.view_dir, cam_param.C);
%         s = R*[start; 1];
%         figure(4);
%         plot3(s(1), s(2), s(3), '*');
%         cam_pos = s(1:3);
%         folder_name = 'rotate_scale_x';
%         applyTransformSaveSVGs(strokes_topology, angle, 1.9, scaleY, scaleZ, focal_point, cam_pos, cam_param, frame_num, folder_name)
%     end
%     
%     frame_num = 0;
%     for i = 0:45:360
%         frame_num = frame_num+1;
%         R=AxelRot(i,cam_param.view_dir, cam_param.C);
%         s = R*[start; 1];
%         figure(4);
%         plot3(s(1), s(2), s(3), '*');
%         cam_pos = s(1:3);
%         folder_name = 'rotate_scale_y';
%         applyTransformSaveSVGs(strokes_topology, angle, scaleX, 1.9, scaleZ, focal_point, cam_pos, cam_param, frame_num, folder_name)
%     end
%     
    
    
%     frame_num = 0;
%     for i = 0:45:360
%         frame_num = frame_num+1;
%         R=AxelRot(i,cam_param.view_dir, cam_param.C);
%         s = R*[start; 1];
%         figure(4);
%         plot3(s(1), s(2), s(3), '*');
%         cam_pos = s(1:3);
%         folder_name = 'back';
%         applyTransformSaveSVGs(strokes_topology, 180, scaleX, scaleY, scaleZ, focal_point, cam_pos, cam_param, frame_num, folder_name)
%     end
    

    function f = distCenter(alpha)
        f = cam_param.C + cam_param.view_dir*alpha';
    end

end





