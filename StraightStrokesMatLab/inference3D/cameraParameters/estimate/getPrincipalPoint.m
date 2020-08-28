function principal_point = getPrincipalPoint(vp, mask_infinite_vps)

global fid;
global DISPLAY_INFO;
%
% Input:
%   vp = 3 by 2 array of image cordinate of vanishign points. [x1 y1; x2
%   y2; x3 y3];
% 
% Output:
%   [x y] cordianted of orthocenter.

if (sum(mask_infinite_vps) == 1)
    % Infiinite vanishign points assumed to be the third:
    principal_point = getPrincipalPoint2fntVP(vp);
    if DISPLAY_INFO
        fprintf(fid, 'Two finit points');
    end
elseif (sum(mask_infinite_vps) == 0)
    principal_point = getPrincipalPoint3fntVP(vp);
    if DISPLAY_INFO
        fprintf(fid, 'Three finit points');
    end
else
    error('Two infinite vanishing points');
end

principal_point = reshape(principal_point, 1, 2);
    
end

function principal_point = getPrincipalPoint2fntVP(vp)
    ref_point = vp(3,:);
    
    r=((ref_point(1)-vp(1,1)).*(vp(2,1)-vp(1,1)) ...
       +(ref_point(2)-vp(1,2)).*(vp(2,2)-vp(1,2)))./...
        ((vp(2,1)-vp(1,1)).^2+(vp(2,2)-vp(1,2)).^2);    
    
    principal_point(1)= vp(1,1) + r.*(vp(2,1)-vp(1,1));
    principal_point(2)= vp(1,2) + r.*(vp(2,2)-vp(1,2));
end

function principal_point = getPrincipalPoint3fntVP(vp)
% Estiamtes prinicpal point assuming three infinite vanishing points.
% The function estimates the thriangle orthocenter, where the triangle is
% given by the three orthogonal vanishing points.
% Finds the slope of two sides of the triangle m
% Finds the slope of two perpenduclars as -1/m
% Express line through opposite vp coordinates and perpendicular slope:
% y - vp_y = (-1/m) * (x - vp_x);
% Solve the system.

% See proof that the prinipal point is an orthocenter:
% "Roberto Cipolla, Tom Drummond, and Duncan P Robertson. 1999. Camera Calibration
% from Vanishing Points in Image of Architectural Scenes.. In BMVC."

side1 = vp(1,:) - vp(2,:);
side2 = vp(1,:) - vp(3,:);

A(1,:) = side1;
A(2,:) = side2;

b(1, 1) = dot(side1, vp(3,:));
b(2, 1) = dot(side2, vp(2,:));

principal_point = A\b;

end
