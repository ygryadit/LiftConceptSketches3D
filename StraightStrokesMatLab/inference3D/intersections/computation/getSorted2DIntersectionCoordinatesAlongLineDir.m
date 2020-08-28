function [p_sorted, ...
          ind_sorted] = getSorted2DIntersectionCoordinatesAlongLineDir(p, pl1, pl2)
    

    num_points = size(p,1);
%     denom = sqrt(sum((pl1  - pl2).^2));
    
    d = (pl2 - pl1)/norm(pl2  - pl1);
    d = repmat(d, num_points, 1);
    v = p-repmat(pl2, num_points, 1);
    t = dot(v,d,2);

    [~,ind_sorted] = sort(t);
    
    p_sorted = p(ind_sorted,:);
  
end