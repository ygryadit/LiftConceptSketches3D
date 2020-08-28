function R = estimateRotationMatrix2VP(vp,principal_point,f)
%  Eq (6) in the paper "Camera calibration using two or three vanishing points" by Orghidan et al.
    principal_point = reshape(principal_point, 1, 2);
    R(:,1) = [vp(1,:)-principal_point f]';
    R(:,1) = R(:,1)/norm(R(:,1));
    R(:,2) = [vp(2,:)-principal_point f]';
    R(:,2) = R(:,2)/norm(R(:,2));
    R(:,3) = cross(R(:,1), R(:,2));    
end
