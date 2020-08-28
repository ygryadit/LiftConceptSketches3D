
function [t, point] = findPointPorjectedPositionOnTheSegment(seg_coordiantes, points_coordiantes)
    
    dx = seg_coordiantes(:,2) - seg_coordiantes(:,1);
    dy = seg_coordiantes(:,4) - seg_coordiantes(:,3);
    
%     lines_lengths = sqrt(dx .^2 + dy .^2);
    
    t = ((points_coordiantes(:,1) - seg_coordiantes(:,1)) .* dx  + ...
         (points_coordiantes(:,2) - seg_coordiantes(:,3)) .* dy) ./ ...
		(dx .* dx + dy .* dy);
    if t < 0
        t = 0;
    end
    
    if t > 1
        t = 1;
    end
    point = [seg_coordiantes(:,1) + t*dx, seg_coordiantes(:,3) + t*dy];
end