function [cam_param, vp, xyz_inds] = ...
                estimateCameraParameters(vp, ...
                                         img_width, ...
                                         point2DConterpartOrigin, ...
                                         img, ...
                                         lambda)
    
    global fid;
    global SHOW_FIGS_PREPROCESS;

    % Axis system:
    % +-1   0   0 -- first vp
    % 0   +-1   0 -- second vp
    % 0     0 +-1 -- third vp

    


    

    [vp_, xyz_inds] = findXYZ(vp);

    mask_infinite_vps = checkInfinityConditionVP(vp_,img_width);
    principal_point_vp = getPrincipalPoint(vp_, mask_infinite_vps);
    
    plotVPsTriangle(1, vp_, principal_point_vp, img);
   

    
    f_vp = estimateFocalLengthFromVPsAndPrincipalPoint(vp_, principal_point_vp);
    fov_vp = focalLength2FieldOfView(f_vp, img_width);
    fprintf(fid, 'Estimated filed of view: fov_vp = %.2f degrees\n', fov_vp);
    
    % Rotation matrix:
    vp = vp_;
    
    if (sum(mask_infinite_vps) == 1)
        R_vp = estimateRotationMatrix2VP(vp, principal_point_vp,f_vp);
        [theta_vp(1), theta_vp(2), theta_vp(3), ~, ~, ~, ~ ] = rotm2eulr(R_vp);
        fprintf(fid, '\t Two finit points\n');
    elseif (sum(mask_infinite_vps) == 0)
        R_vp = estimateRotationMatrix3VP(vp,f_vp, principal_point_vp);
        fprintf(fid, '\t Three finit points\n');
    else
        fprintf(fid, '\t Two vanishing points at infinity: near orthogonal projection\n');
        R_vp = eye(3)*NaN;
%         C_vp = ones(3,1)*NaN;
        return;
    end
    
    if ~exist('lambda', 'var')
        lambda = 1.0;
    end
    
    t_vp = translationUpToScale(lambda, principal_point_vp, ...
                                      point2DConterpartOrigin, ...
                                      f_vp);
     t_vp = reshape(t_vp, 3, 1);                             
%     t_vp = computeTranslationVector(principal_point_vp, ...
%                                     point2DConterpartOrigin, ...
%                                       R_vp, f_vp,...
%                                       point2D, point3D);
    
    % Calibration matrix:
      K_vp = [f_vp   0       principal_point_vp(1) 0; 
             0     f_vp    principal_point_vp(2) 0;
             0     0       1.0                   0];
         
    T_vp = [R_vp t_vp; zeros(1,3) 1];
    P_vp = K_vp*T_vp;
    
    
    p1 = P_vp*[1 0 0 1]';
    p1 = p1/p1(3);
    p2 = P_vp*[0 1 0 1]';
        p2 = p2/p2(3);
    p3 = P_vp*[0 0 1 1]';
        p3 = p3/p3(3);
    
        
    if SHOW_FIGS_PREPROCESS
        
        plot([point2DConterpartOrigin(1) p1(1)], [point2DConterpartOrigin(2) p1(2)],'m-*');
        plot([point2DConterpartOrigin(1) p2(1)], [point2DConterpartOrigin(2) p2(2)],'g-*');
        plot([point2DConterpartOrigin(1) p3(1)], [point2DConterpartOrigin(2) p3(2)],'b-*');
        text(p1(1), p1(2),'x');
        text(p2(1), p2(2),'y');
        text(p3(1), p3(2),'z');
        axis equal;
        axis ij;
        
    end
    
    
  
    
    
    cam_param.P = P_vp;
    cam_param.f = f_vp;
    cam_param.principal_point = principal_point_vp;
    cam_param.R = R_vp;
    cam_param.K = K_vp(1:3, 1:3);
    cam_param.t = t_vp;
    cam_param.view_dir = T_vp(3, 1:3);
    cam_param.fov = fov_vp;
    cam_param.C = -inv(R_vp)*t_vp;
%     p1 = P_vp*[0 0 0 1]';
%     p1 = p1/p1(3);
%     disp(p1(1:2));
%     disp(point2DConterpartOrigin);
%     
%     p2 = P_vp*[point3D 1]';
%     p2 = p2/p2(3);
%     disp(p2(1:2));
%     disp(point2D);


end

function [vp_, xyz_inds] = findXYZ(vp)
    %Vanishing points on a new order:
    vp_ = zeros(3,2);
    xyz_inds = zeros(3,1);
    
    
    % z vanishing point:
    [~,vals] = sort(vp(:,1));
    z_pos = vals(2);
   
    %Vertical direction
    vp_(3, :) = vp(z_pos, :);  
    xyz_inds(3) = z_pos;
    
    
    % x,y vanishing points:
    inds = vals([1,3]);
    [~, left_vp] = max(vp_(3, 1) - vp(inds, 1));
    [~, right_vp] = min(vp_(3, 1) - vp(inds, 1));
    
    left_vp = inds(left_vp);
    right_vp = inds(right_vp);
    mean_y = mean(vp([left_vp,right_vp],2));
    
       
    global ZUP;
    ZUP = vp(z_pos,2) < mean_y;   
    
    %Rigth coordinate system:
    if ZUP
        %z is up
        %left vp is 'y' and right is 'x'
        vp_(2, :) = vp(left_vp, :);
        xyz_inds(2) = left_vp;
        vp_(1, :) = vp(right_vp, :);
        xyz_inds(1) = right_vp;
    else
        %z is down
        %left vp is 'x' and right is 'y'
        vp_(1, :) = vp(left_vp, :);
        xyz_inds(1) = left_vp;
        vp_(2, :) = vp(right_vp, :);
        xyz_inds(2) = right_vp;
    end
        
end


function plotVPsTriangle(fig_num, vp, principal_point_vp, img)
    global SHOW_FIGS_PREPROCESS;
    if SHOW_FIGS_PREPROCESS
        close(figure(fig_num));
        figure(fig_num);
        plot(vp(:,1), vp(:,2),'*');
        hold on;
        imagesc(img);
        axis equal;
        axis ij;

        %Vanishig points:
        plot(vp(1,1), vp(1,2),'*');
        plot(vp(2,1), vp(2,2),'*');
        plot(vp(3,1), vp(3,2),'*');
        text(vp(1,1), vp(1,2),'x');
        text(vp(2,1), vp(2,2),'y');
        text(vp(3,1), vp(3,2),'z');
        %Vanishing lines:
        plot(vp(1:2,1), vp(1:2,2));
        plot(vp(1:3,1), vp(1:3,2));
        plot(vp(2:3,1), vp(2:3,2));

        %Principle point:
        plot(principal_point_vp(:,1), principal_point_vp(:,2),'o');
        
    end
    

end