function orthoViews(strokes_topology, intersections, cam_param)
    focal_point = [0,0,0];
    up = [0,0,1];
    cam_pos = [1, 0, 0];   
    cam_pos = reshape(cam_pos, 1, 3);
    cam_param.R = rotationMatrixFromView(cam_pos, focal_point, up);
    cam_param.C = cam_pos';
    cam_param.K(1,3) = 512/2;
    cam_param.K(2,3) = 512/2;

    P = cam_param.K*[ cam_param.R -cam_param.R*cam_param.C];
    
        
    img = ones(512,512,3);
    cam_param.P = P;    
    
    reproject3Dto2D(img, cam_param, strokes_topology,[],2);
end