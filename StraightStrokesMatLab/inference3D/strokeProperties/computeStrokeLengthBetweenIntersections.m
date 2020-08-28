function length2D = computeStrokeLengthBetweenIntersections(stroke, ind_inter1, ind_inter2)


ind_seg1 = stroke.inter_seg_nums(ind_inter1);
ind_seg2 = stroke.inter_seg_nums(ind_inter2);

[point,v1,v2] = getPoint(ind_seg1, stroke.poly2d_extended);
dist1 = norm(v2 - point);

[point,v1,v2] = getPoint(ind_seg2, stroke.poly2d_extended);
dist2 = norm(point - v1);
    
l = getLengthPolyVertexVertex(stroke.poly2d_extended, ...
                              ceil(ind_seg1), ...
                              floor(ind_seg2));

length2D = dist1 + dist2 + l;

end

function [point,v1,v2] = getPoint(ind_seg, poly)
ind_vertex_1 = floor(ind_seg);
ind_vertex_2 = ceil(ind_seg);
seg_pos = ind_seg - ind_vertex_1;

v1 = poly(ind_vertex_1,:);
v2 = poly(ind_vertex_2,:);

delta = v2 - v1;
point = v1 + delta*seg_pos;
end

function l = getLengthPolyVertexVertex(poly, vertex1, vertex2)
poly = poly(vertex1:vertex2,:);
l = 0;
if size(poly,1) == 0
    return;
end
for i = 1:(size(poly,1)-1)
  l = l + norm(poly(i+1,:) - poly(i,:)); 
end
end