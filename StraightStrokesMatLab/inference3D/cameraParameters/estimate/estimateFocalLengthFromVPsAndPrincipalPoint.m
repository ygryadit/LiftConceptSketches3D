function f = estimateFocalLengthFromVPsAndPrincipalPoint(vp, principal_point)
% "Camera calibration using two or three vanishing points" by Orghidan et
% al. Eq.(18)
% Input:
%   vp: 3x2 --- coordiantes of three vanishing points:
principal_point = reshape(principal_point, [1,2]);


if (size(vp, 2) > 3)
    vp = reshape(vp, 2,3)';
end


f = sqrt(abs(dot( vp(1, :) - principal_point, principal_point - vp(2, :))));

end