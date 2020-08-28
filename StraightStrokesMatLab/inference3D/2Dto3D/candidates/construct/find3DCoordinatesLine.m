function coordiantes3d = find3DCoordinatesLine(P, line_origin, line_dir, coordiantes2d)
u = coordiantes2d(1);
v = coordiantes2d(2);
A = [dot(P(1,1:3) - u*P(3,1:3), line_dir);...
     dot(P(2,1:3) - v*P(3,1:3), line_dir)];   
b = [u*P(3,4) - P(1,4) + dot( u*P(3,1:3) - P(1,1:3),  line_origin);... 
     v*P(3,4) - P(2,4) + dot( v*P(3,1:3) - P(2,1:3),  line_origin)];

 t = A\b;
 coordiantes3d = line_origin + t*line_dir;
end