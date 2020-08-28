% function [distances, lines_lengths] = findLinePointDistance(lines_coordinates, points_coordiantes)
% 
% Finds the distane beetween a point and a line segment
% 
% Old: Returns the coordinate of the points in between two strokesm this value should be recomputed precisely looking at polylines.
% New: Returns the coordinate of the two points on the two strokes this value should be recomputed precisely looking at polylines.

function [distances, intersection_coordinates] = find2DLineSegmentPointDistance(lines_coordinates, points_coordiantes)
    
    dx = lines_coordinates(:,2) - lines_coordinates(:,1);
    dy = lines_coordinates(:,4) - lines_coordinates(:,3);
    
%     lines_lengths = sqrt(dx .^2 + dy .^2);
    
    t = ((points_coordiantes(:,1) - lines_coordinates(:,1)) .* dx  + ...
         (points_coordiantes(:,2) - lines_coordinates(:,3)) .* dy) ./ ...
		(dx .* dx + dy .* dy);
    
    ind_end1 = find(t < 0);
    ind_end2 = find(t > 1);
    ind_line = find(t >= 0 & t <=1);
    
    distances(ind_end1) = sqrt( (points_coordiantes(ind_end1,1) - lines_coordinates(ind_end1,1)).^2 + ...
                                (points_coordiantes(ind_end1,2) - lines_coordinates(ind_end1,3)).^2 );
    distances(ind_end2) = sqrt( (points_coordiantes(ind_end2,1) - lines_coordinates(ind_end2,2)).^2 + ...
                                (points_coordiantes(ind_end2,2) - lines_coordinates(ind_end2,4)).^2 );
    distances(ind_line) = sqrt( (points_coordiantes(ind_line,1) - (lines_coordinates(ind_line,1) + t(ind_line).* dx(ind_line)  ) ).^2 + ...
                                (points_coordiantes(ind_line,2) - (lines_coordinates(ind_line,3) + t(ind_line).* dy(ind_line)  ) ).^2 );
    distances = reshape(distances, [], 1);                   
    
    
    intersection_coordinates = zeros(size(lines_coordinates,1), 4);

%     intersection_coordinates(ind_end1, :) = 0.5*(points_coordiantes(ind_end1,:) + lines_coordinates(ind_end1,[1,3]));
%     intersection_coordinates(ind_end2, :) = 0.5*(points_coordiantes(ind_end2,:) + lines_coordinates(ind_end2,[2,4]));
%     intersection_coordinates(ind_line, 1) = 0.5*(points_coordiantes(ind_line,1) + (lines_coordinates(ind_line,1) + t(ind_line).* dx(ind_line)));
%     intersection_coordinates(ind_line, 2) = 0.5*(points_coordiantes(ind_line,2) + (lines_coordinates(ind_line,3) + t(ind_line).* dy(ind_line)));
    

    intersection_coordinates(ind_end1, :) = [points_coordiantes(ind_end1,:)  lines_coordinates(ind_end1,[1,3])];
    intersection_coordinates(ind_end2, :) = [points_coordiantes(ind_end2,:)  lines_coordinates(ind_end2,[2,4])];
    intersection_coordinates(ind_line, [1 3]) = [points_coordiantes(ind_line,1) (lines_coordinates(ind_line,1) + t(ind_line).* dx(ind_line))];
    intersection_coordinates(ind_line, [2 4]) = [points_coordiantes(ind_line,2) (lines_coordinates(ind_line,3) + t(ind_line).* dy(ind_line))];



%     intersection_coordinates(ind_line, :) = 0.5*(points_coordiantes(ind_line,:) + ...
%                                                 (lines_coordinates(ind_line,[1,3]) + repmat(t(ind_line).* dy(ind_line), 1,2)  ) ...
%                                                  );
    

    
end