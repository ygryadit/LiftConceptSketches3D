function view_matrix = rotationMatrixFromView(cam_pos, focal_point, up)
    view_matrix = eye(3);
    direction = ( focal_point - cam_pos);
    direction = direction./norm(direction);
    
    s = cross(-direction, up);
    s = s/norm(s);
    
    u = cross(s, -direction);
    u = u/norm(u);
    
   
    view_matrix(1,1:3) = s;
    view_matrix(2,1:3) = u;
    view_matrix(3,1:3) = direction;
%     view_matrix(1,4) = -s*cam_pos';
%     view_matrix(2,4) = -u*cam_pos';
%     view_matrix(3,4) = direction*cam_pos';
end