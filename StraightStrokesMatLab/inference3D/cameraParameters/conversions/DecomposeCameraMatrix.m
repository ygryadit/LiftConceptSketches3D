function [k,r,c] = DecomposeCameraMatrix(p) 
% Descompose the camera matrix P into K, R, C according P = K [R |-RC]
%
% Input:
%       - p is a 3x4 camera matrix.
%
% Output:
%       - k is a 3x3 calibration matrix.
%       - r is a 3x3 rotation matrix.
%       - c is a 3x1 camera coordinates.
%
%----------------------------------------------------------
%
% From 
%    Book: "Multiple View Geometry in Computer Vision",
% Authors: Hartley and Zisserman, 2006, [HaZ2006]
% Section: "Decomposition of the camera matrix", 
% Chapter: 6
%    Page: 163
%
%----------------------------------------------------------
%      Author: Diego Cheda
% Affiliation: CVC - UAB
%        Date: 03/06/2008
%----------------------------------------------------------

% Finding camera centre
% PC = 0, right null vector of P.
x =  det([ p(:,2), p(:,3), p(:,4) ]);
y = -det([ p(:,1), p(:,3), p(:,4) ]);
z =  det([ p(:,1), p(:,2), p(:,4) ]);
t = -det([ p(:,1), p(:,2), p(:,3) ]);

c = [ x/t; y/t; z/t ];

% Finding camera orientation and internal parameters
% P = [M | -MC] = K[R | -RC]
[k,r] = RQDecomposition(p(:,1:3));


% for (i=1:3)
%    for (j=1:3)
%         if (i==j) && (k(i,j) < -1e-5)
%             t = k;
%             k(i,j) = k(i,j)*-1;
%             m = t \ k;
%             r = m\r;
%         end
%    end
% end

% disp(k*r)