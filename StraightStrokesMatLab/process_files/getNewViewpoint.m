function getNewViewpoint()
%     folder1 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400\Prof2task2\guitar_01\view1';
%     folder2 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_all_inter_400\_400\Prof2task2\guitar_01\view1';

    folder1 = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\curves17052020\results\results_good';
    folder2 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400_no_delay\_400\Professional2\vacuum_cleaner\view1';
    
%     filename = 'Prof2task2_guitar_01';
    filename = 'Professional2_vacuum_cleaner';
    
    global folder_save;
    folder_save = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\views';
    if ~exist(folder_save, 'dir')
        mkdir(folder_save);
    end
%       [ strokes_topology, intersections, cam_param] = ...
%             readReconstructionJsonFelix( fullfile(folder, filename)); 
 
    [ strokes_topology1, intersections, cam_param] = ...
        readReconstructionJsonFelix( fullfile(folder1, [filename '_bestScore_full.json'] ));

    [strokes_topology1, intersections, cam_param] = ...
        scaleAndCenterTheObject(strokes_topology1, intersections, cam_param, []);   

    [ strokes_topology2, intersections, cam_param] = ...
        readReconstructionJson( fullfile(folder2, [filename '_final_full.json'] ));

    [strokes_topology2, intersections, cam_param] = ...
        scaleAndCenterTheObject(strokes_topology2, intersections, cam_param, []);   
    
    angle = 0;
    scaleX = 1.0;
    scaleY = 1.0;
    scaleZ = 1.0;
    
    cam_param.view_dir = cam_param.view_dir./norm(cam_param.view_dir);

    alpha = lsqnonlin(@distCenter, norm(cam_param.view_dir+cam_param.C));
    focal_point = cam_param.C + alpha * cam_param.view_dir;

    frame_num = 0;
    folder = filename;
    
    folder_name = [folder '_v49planes_400'];    
    applyTransformSaveSVGs(strokes_topology1, angle, 0.7, scaleY, scaleZ, focal_point, cam_param.C, cam_param, frame_num, folder_name, intersections)
    
%     folder_name = [folder '_v49planes_all_inter_400'];    
%     applyTransformSaveSVGs(strokes_topology2, angle, scaleX, scaleY, scaleZ, focal_point, cam_param.C, cam_param, frame_num, folder_name, intersections)
%  
    frame_num = 1;
    folder_name = [folder '_v49planes_400'];    
    applyTransformSaveSVGs(strokes_topology1, angle, scaleX, 0.5, scaleZ, focal_point, cam_param.C, cam_param, frame_num, folder_name, intersections)
    
%     folder_name = [folder '_v49planes_all_inter_400'];    
%     applyTransformSaveSVGs(strokes_topology2, angle, scaleX, scaleY, scaleZ, focal_point, cam_param.C, cam_param, frame_num, folder_name, intersections)
    
    frame_num = 2;
    cam_pos = cam_param.C - [0, 0, 0.0]';
    folder_name = [folder '_v49planes_400'];    
    applyTransformSaveSVGs(strokes_topology1, angle, 0.7, scaleY, scaleZ, focal_point, cam_pos , cam_param, frame_num, folder_name, intersections)
    
%     folder_name = [folder '_v49planes_all_inter_400'];    
%     applyTransformSaveSVGs(strokes_topology2, angle, scaleX, scaleY, scaleZ, focal_point, cam_pos , cam_param, frame_num, folder_name, intersections)
    
    
    function f = distCenter(alpha)
        f = cam_param.C + cam_param.view_dir*alpha';
    end

end