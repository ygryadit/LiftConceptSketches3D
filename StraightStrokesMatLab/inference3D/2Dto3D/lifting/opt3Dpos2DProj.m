function [points3D, lambdas] = opt3Dpos2DProj(points_2D, cam_param, ...
                                              p_line_prior_3D_1, p_line_prior_3D_2)
% Optimise to 3D position according to a 2D projection and 3D prior:
% 
% Input:
%   points2D : N x 2, N - number of points
% 

    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
                            'Display', 'off');  
    
    num_points = size(points_2D,1);
    lambdas = ones(num_points, 1);
    lb =[];
    ub =[];
    
    [lambdas,...
     resnorm,...
     residual,...
     exitflag,...
     output] = lsqnonlin(@distPointGeomPrior, lambdas, lb, ub, options);
    
    points3D = lambdasTo3DCoordinates(lambdas, points_2D, cam_param);
    
    %% Function used in the optimisation:
    function distances = distPointGeomPrior(lambdas)
        points_3D = lambdasTo3DCoordinates(lambdas, points_2D, cam_param);
    
        distances = pointLineDistance3D(points_3D, p_line_prior_3D_1, p_line_prior_3D_2)';
    end

end