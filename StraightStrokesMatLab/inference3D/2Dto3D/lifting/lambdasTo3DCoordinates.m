% function bezier_coordinates3D = lambdasTo3DCoordinates(bezier_lambdas, bezier_coordinates_2D, cam_param)
% 
% Takes as input projection matrix, 2D coordinate of the point and a scale
% factor \lamda. Returns the corresponidng 3D coordinate.

function bezier_coordinates3D = lambdasTo3DCoordinates(bezier_lambdas, bezier_coordinates_2D, cam_param)

    num_points = length(bezier_lambdas);
    bezier_coordinates3D = zeros(num_points,3);
    R = cam_param.R;
    t = cam_param.t;
    f = cam_param.f;
    u0 = cam_param.principal_point(1);
    v0 = cam_param.principal_point(2);
      
    K = [f 0 u0; 0 f v0; 0 0 1];
    for i = 1:num_points
        lambda = bezier_lambdas(i);
        u = bezier_coordinates_2D(i,1);
        v = bezier_coordinates_2D(i,2);
        
        p_cam = inv(K)*[u; v; 1];
        
        p_world = inv([R t; 0 0 0 1])*[p_cam*lambda; 1];
        bezier_coordinates3D(i,:) = p_world(1:3);
        
%         vec = cam_param.P*[bezier_coordinates3D(i,:)'; 1];   
%         vec(1:2) = vec(1:2)/vec(3);
%         fprintf('(u,v) = (%.2f, %.2f), estiamted ((%.2f, %.2f))\n', u,v, vec(1:2));
    end
% 
%     bezier_coordinates3D(:,1) = bezier_lambdas.*(bezier_coordinates_2D(:,1) - cam_param.principal_point(1))/cam_param.ff;
%     bezier_coordinates3D(:,2) = bezier_lambdas.*(bezier_coordinates_2D(:,2) - cam_param.principal_point(2))/cam_param.ff;
%     bezier_coordinates3D(:,3) = bezier_lambdas;
end
