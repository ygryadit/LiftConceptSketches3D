function q = GivensRotationMatrix(matrixA, axis)
% Computes the Givens rotation matrix.
% 
% Input:
%       - matrixA is matrixA 3x3 matrix.
%
%       - axis is the rotation over one of the 3 coordinate axes (1 for x, 2
%       for y, and 3 for z).
%
% Output:
%       - q is the rotation matrix.
%
%----------------------------------------------------------
%
% From 
%    Book: "Multiple View Geometry in Computer Vision",
% Authors: Hartley and Zisserman, 2006, [HaZ2006]
% Section: "Givens rotations and RQ decomposition", 
% Chapter: appendix 4
%    Page: 579
%
%----------------------------------------------------------
%      Author: Diego Cheda
% Affiliation: CVC - UAB
%        Date: 03/06/2008
%----------------------------------------------------------

if (axis == 1)
    % over x
	%  _        _
	% | 1        |
	% |    c  -s |
	% |    s   c |
	%  -        -
	c = -matrixA(3,3)/sqrt(matrixA(3,2)^2 + matrixA(3,3)^2);
	s = matrixA(3,2)/sqrt(matrixA(3,2)^2 + matrixA(3,3)^2);
	q = [1 0 0; 0 c -s; 0 s c];
elseif (axis == 2)
	% over y
	%  _        _
	% | c      s |
	% |    1     |
	% | -s     c |
	%  -        -
	c = matrixA(3,3)/sqrt(matrixA(3,1)^2 + matrixA(3,3)^2);
	s = matrixA(3,1)/sqrt(matrixA(3,1)^2 + matrixA(3,3)^2);
	q = [c 0 s; 0 1 0; -s 0 c];
elseif (axis == 3)
	% over z
	%  _        _
	% | c -s     |
	% | s  c     |
	% |        1 |
	%  -        -
	c = -matrixA(2,2)/sqrt(matrixA(2,2)^2 + matrixA(2,1)^2);
	s = matrixA(2,1)/sqrt(matrixA(2,2)^2 + matrixA(2,1)^2);
	q = [c -s 0; s c 0; 0 0 1];
end
