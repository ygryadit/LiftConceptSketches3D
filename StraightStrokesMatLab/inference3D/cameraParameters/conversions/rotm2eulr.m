function [ x1,y1,z1,  x2,y2,z2, R] = rotm2eulr( R )
% R = R_z*R_y*R_x;
% 
% [ cos(theta_y)*cos(theta_z), cos(theta_z)*sin(theta_x)*sin(theta_y) - cos(theta_x)*sin(theta_z), sin(theta_x)*sin(theta_z) + cos(theta_x)*cos(theta_z)*sin(theta_y)]
% [ cos(theta_y)*sin(theta_z), cos(theta_x)*cos(theta_z) + sin(theta_x)*sin(theta_y)*sin(theta_z), cos(theta_x)*sin(theta_y)*sin(theta_z) - cos(theta_z)*sin(theta_x)]
% [             -sin(theta_y),                                          cos(theta_y)*sin(theta_x),                                          cos(theta_y)*cos(theta_x)]
% 
% theta_y=asin(R_est(1,3));
% theta_x=acos(min( max(R_est(1,1)/cos(theta_y), -1) , 1.0));
% 
% theta_z=acos(min( max(R_est(3,3)/cos(theta_y), -1) , 1.0));
% 
% rx = theta_x/pi*180
% ry = theta_y/pi*180
% rz = theta_z/pi*180

% sy = sqrt(R(1,1).^2 +  R(2,1).^2);
%  
% singular = sy < 1e-6; 
%  
%     
% if ~singular
% %     x = atan2(R(3,2), R(3,3));
% %     if (cos(x) == 0)
% %         y = atan2(-R(3,1), R(3,2)/sin(x));
% %     else
% %         y = atan2(-R(3,1), R(3,3)/cos(x));
% %     end
% %     z = atan2(R(2,1), R(1,1));
%     
% %     y1 = asin(-R(3,1));
%     y1 = -asin(R(3,1));
%     y2 = pi - y1;
%     
%     x1 = atan2(R(3,2)/cos(y1), R(3,3)/cos(y1));
%     x2 = atan2(R(3,2)/cos(y2), R(3,3)/cos(y2));
%     
%     z1 = atan2(R(2,1)/cos(y1), R(1,1)/cos(y1));
%     z2 = atan2(R(2,1)/cos(y2), R(1,1)/cos(y2));
%     
%     
%     
% else
% % [ 0, cos(theta_z)*sin(theta_x)*sin(theta_y) - cos(theta_x)*sin(theta_z), sin(theta_x)*sin(theta_z) + cos(theta_x)*cos(theta_z)*sin(theta_y)]
% % [ 0, cos(theta_z)*cos(theta_x) + sin(theta_x)*sin(theta_y)*sin(theta_z), cos(theta_x)*sin(theta_y)*sin(theta_z) - cos(theta_z)*sin(theta_x)]
% % [             -sin(theta_y),                                          0,                                          0]
% %     x = atan2(-R(2,3), R(2,2));
% %     y = atan2(-R(3,1), sy);
% %     z = 0;
%     y = asin(-R_est(3,1));
%     z = 0;
%     if (R_est(3,1) < 0)
%        x = y + atan2(R(1,2), R(1,3));
%     else
%        x = -y + atan2(-R(1,2), -R(1,3));
%     end
% end




 
    
if (abs(R(3,1)) ~= 1)

    y1 = -asin(R(3,1));
    y2 = pi - y1;
    
    x1 = atan2(R(3,2)/cos(y1), R(3,3)/cos(y1));
    x2 = atan2(R(3,2)/cos(y2), R(3,3)/cos(y2));
    
    z1 = atan2(R(2,1)/cos(y1), R(1,1)/cos(y1));
    z2 = atan2(R(2,1)/cos(y2), R(1,1)/cos(y2));
    
else
    z1 = 0;
    z2 = 0;
    
    if (R(3,1) == -1)
        y1 = pi/2;
        y2 = y1;
        
        x1 = z1 +  atan2(R(1,2), R(1,3));
        x2 = x1;
%        x = y + atan2(R(1,2), R(1,3));
    else
        y1 = -pi/2;
        y2 = y1;
        x1 = -z1 + atan2(-R(1,2), -R(1,3));
        x2 = x1;
    end
end

theta_z = z1;
theta_x = x1;
theta_y = y1;

    R_z =    [cos(theta_z) -sin(theta_z) 0;
           sin(theta_z)  cos(theta_z) 0;
            0 0 1];

    R_y = [cos(theta_y)   0 sin(theta_y);
                0            1 0;
           -sin(theta_y)  0 cos(theta_y)];


    R_x = [ 1 0 0; 
            0 cos(theta_x) -sin(theta_x);
            0 sin(theta_x)  cos(theta_x)];


    R =  R_z*R_y*R_x;     

end

