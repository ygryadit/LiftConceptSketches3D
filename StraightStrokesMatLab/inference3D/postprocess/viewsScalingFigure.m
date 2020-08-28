function viewsScalingFigure()

folder = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\json_cluster_v3';
global folder_save;
global sketch_height;
sketch_height = 512;

files = dir(folder);
files = files(~[files(:).isdir]);
files = {files.name};


global sketch_height;
sketch_height = 512;

folder = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\json_cluster_v3';
foler_save_base = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\scaling_figure\';

for i = [1 5]
   filename = files{i};
   filename_ = strrep(filename, '_bestScore_full.json', '');
% viewsObject('designer2_guitar_01', foler_save_base, folder)
% viewsObject('Prof2task2_cabinet_01', foler_save_base, folder)
    close all
    viewsObject(filename_, foler_save_base, folder)
end

end





function viewsObject(filename, foler_save_base, folder)
    global folder_save;
    try
    [ strokes_topology, intersections, cam_param] = ...
        readReconstructionJsonFelix( fullfile(folder, [filename '_bestScore_full.json'] ));
    catch
        return;
    end
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
    scaleY = 2.0;
    scaleZ = 1.0;
    
    
     folder_name = 'back';
   applyTransformSaveSVGs(strokes_topology, angle, scaleX, scaleY, scaleZ, focal_point, cam_param.C, cam_param, frame_num, folder_name)

   
    
    
    
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





