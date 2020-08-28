function distances = pointLineDistance3D(p, pl1, pl2)
    num_points = size(p,1);
%     denom = sqrt(sum((pl1  - pl2).^2));
    
    d = (pl2 - pl1)/norm(pl2  - pl1);
    d = repmat(d, num_points, 1);
    v = p-repmat(pl2, num_points, 1);
    t = dot(v,d,2);
    P = repmat(pl2, num_points, 1) + t.*d;
    
    distances = P - p;
    distances = sqrt(sum(distances.^2,2));
end