% Carsten Rother: "A new approach to vanishing point detection in architectural
% environments"
% Section 3.1:
% Computes the distance between the line segment 'segment_s' and the line
% passing through the midpoint of the 'segment_s' and vanishing points 'vp'.
% 
% Output:
%   theta: an angle in radians between the line segment 'segment_s' and the 
%       line passing through the midpoint of the 'segment_s' and vanishing 
%       points 'vp'.
%   vp_outside_line_segment: boolean, to insure that the line point is  not
%       liying withing line segment endpoints.
%       Since in sketches sometimes actuall lines towards vanishing points
%       are drawn the condition is modified to 0.95
function   [theta, vp_outside_line_segment]=computeDistanceLineVPs(segment_s, vp)
% Line segment s (Fig 3 a)
x1=segment_s(1); 
x2=segment_s(2);
y1=segment_s(3);
y2=segment_s(4);

midpnt=[(x1+x2)/2  (y1+y2)/2];
lengthl=sqrt((x1-x2)^2+(y1-y2)^2);

slope1=(y2-y1)/(x2-x1);

% Line segment s_dash (Fig 3 a)
slope2=(vp(:,2)-ones(size(vp,1),1)*midpnt(2))./...
       (vp(:,1)-ones(size(vp,1),1)*midpnt(1));

% Angle between s and s_dash
slope1 = ones(size(slope2,1),1)*slope1;
theta=atan(abs((slope1-slope2(:))./(1+slope1.*slope2(:))));

%Discard vanishing points that lie withing the endpoints of s' (segment_s
%rotated to lie on line towards vanishing point)
% Distance between midpoint and vanishing point:
d=sqrt((vp(:,1)-ones(size(vp,1),1)*midpnt(1)).^2+(vp(:,2)-ones(size(vp,1),1)*midpnt(2)).^2);
vp_outside_line_segment = false(size(vp,1),1);
vp_outside_line_segment(d >= 0.95*lengthl/2) = true;
return

