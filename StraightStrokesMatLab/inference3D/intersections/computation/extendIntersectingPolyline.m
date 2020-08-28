function [begin_x, begin_y, end_x, end_y] = extendIntersectingPolyline(sketch_stroke)
        dir(1)= sketch_stroke.primitive_geom(:,2) - sketch_stroke.primitive_geom(:,1);
        dir(2) = sketch_stroke.primitive_geom(:,4) - sketch_stroke.primitive_geom(:,3);
    
%         dir(1) = 
%         dir_(1) = sketch_strokes(s1).points2D(2).x - sketch_strokes(s1).points2D(1).x;
%         dir(2) = sketch_strokes(s1).points2D(2).y - sketch_strokes(s1).points2D(1).y;
        dir = dir/norm(dir);
        
        begin_x = sketch_stroke.primitive_geom(:,1)  - sketch_stroke.accuracy_radius*dir(1);
        begin_y = sketch_stroke.primitive_geom(:,3)  - sketch_stroke.accuracy_radius*dir(2);
        
        end_x = sketch_stroke.primitive_geom(:,2)  + sketch_stroke.accuracy_radius*dir(1);
        end_y = sketch_stroke.primitive_geom(:,4)  + sketch_stroke.accuracy_radius*dir(2);      
end