function pAngle = ...
                probabilitiesAngle( lines_dirs, direction_prior)
% Description:
% Map the angles so that zero difference gives probablility one and 
% 0.14 = 1 - cos(pi/6) - the value where the weight should be equal to zero.

global sigma_aligned;

try
    lines_dirs = normilizeDirections(lines_dirs);
catch e
    rethrow(e);
end
    
num_lines = size(lines_dirs,1);
direction_prior = repmat(direction_prior, num_lines, 1);

cos_angle = abs(dot(lines_dirs, direction_prior, 2));


% 0.0044 = 2*(0.14/3)^2
% sigma_old = 0.05;
%sigma = sin(pi/36)/3;
% sigma = sin(pi/45);
% sigma = sind(1);
%sigma = sind(0.75);

% sigma = sin(pi/180)/3;
% sigma = sin(pi/90)/3;
pAngle = exp(-(1.0 - cos_angle).^2/(2*(sigma_aligned)^2));
end


function lines_dirs = normilizeDirections(lines_dirs)
    lengths     = repmat(sqrt(sum(lines_dirs.^2, 2)), 1, 3);
    lines_dirs  = lines_dirs./lengths ;
end






