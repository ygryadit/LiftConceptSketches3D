function applyTransformSaveSVGs(strokes_topology, angle, scaleX, scaleY, scaleZ, focal_point, cam_pos, cam_param, frame_num, file_name, intersections)

    global save_color;
    % Folder
%     [folder_black, folder_color] = setupFOlder(folder_name);
    global folder_save;
    
    folder_black = folder_save;
    folder_color = folder_save;

    if ~exist(folder_black, 'dir')
       mkdir(folder_black);
    end
    
    % File paths
    file_svg = fullfile(folder_black, sprintf('%s%03d.svg', file_name, frame_num))
    if save_color 
        file_svg_color = fullfile(folder_color, sprintf('sketch_color%03d.svg', frame_num));
    else
        file_svg_color = '';
    end
    
   
    % Update camera matrix:
    cam_param.P
    P = updateCameraMatrix(cam_param, focal_point, cam_pos)

    % Object rotation:
    R_o = rotz_mine(angle);
    
    % Project and save
    strokes_topology = projectPointsSaveSVG(strokes_topology, scaleX, scaleY, scaleZ, R_o, P, file_svg, file_svg_color);
    
    img = ones(512,512,3);
    cam_param.P = P;    
%     reproject3Dto2D(img, cam_param, strokes_topology,[],2,'k');
    
    
    object.cam_param    = cam_param;
    object.strokes_topology = strokes_topology;
    object.intersections = intersections;
    object_str = jsonencode(object);

    file_json = fullfile(folder_black, sprintf('%s%03d.json', file_name, frame_num));
        
    fid_ = fopen(file_json, 'w');
    fwrite(fid_, object_str);
    fclose(fid_); 
end

function P = updateCameraMatrix(cam_param, focal_point, cam_pos)
    up = [0 0 1];
    cam_param.C = cam_pos;
    cam_pos = reshape(cam_pos, 1, 3);
    cam_param.R = rotationMatrixFromView(cam_pos, focal_point', cam_param.R(2,:));
    P =  cam_param.K *[ cam_param.R -cam_param.R*cam_param.C];
end
