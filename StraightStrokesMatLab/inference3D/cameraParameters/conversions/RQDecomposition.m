function [k,r] = RQDecomposition(matrixA)
% Computes the RQ decomposition.
%
% Input:
%       - a is a 3x3 matrix.
%
% Output:
%       - k is a 3x3 upper triangular matrix.
%           
%       - r is a 3x3 rotation matrix.
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


% In
% addition, the angles ?x, ?y and ?z associated with the three Givens rotations provide a
% parametrization of the rotation by three Euler angles, otherwise known as roll, pitch
% and yaw angles

    qx = GivensRotationMatrix(matrixA,1);

    k = matrixA * qx;

    qy = GivensRotationMatrix(k,2);

    k = k * qy;

    qz = GivensRotationMatrix(k,3);

    k = k * qz;

    r = qz' * qy' * qx';

    if ((k(1,1) < 0) && (k(2,2) < 0))
        d = diag([-1 -1 ones(1,1)]);
        k = k * d;
        r = d * r;
    end

%     if k(3,3) > 0
%         d = diag([1 1 -1]);
%         k = k * d;
%         r = d * r;
%     end

%     if k(1,1) < 0
%         
%     end
end