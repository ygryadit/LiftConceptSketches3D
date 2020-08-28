function R = estimateRotationMatrix3VP(vp,f, principal_point)
%   Recover scaling parameters of vanishing points in Eq. (14) using Eq. (22)
%   Implements paper "Camera calibration using two or three vanishing points" by Orghidan et al.
% 	Section III. Camera calibratino using three vanishing points.
    A(1,:) = vp(:,1)';
    A(2,:) = vp(:,2)';
    A(3,:) = vp(:,1).^2';
    A(4,:) = vp(:,2).^2';
    A(5,:) = vp(:,1)'.*vp(:,2)';
    
    b(1:2) = principal_point;
    b(3:4)= f^2 + principal_point.^2;
    b(5) = principal_point(1)*principal_point(2);
    
    lambda = sqrt(A\b');
   
    R(1,:) = lambda.*(vp(:,1) - principal_point(1))/f;
    R(2,:) = lambda.*(vp(:,2) - principal_point(2))/f;
    R(3,:) = lambda;
end
