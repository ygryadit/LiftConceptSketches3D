% -------------------------------------------------------------------------
% Input:
% -------------------------------------------------------------------------
% ------------------------------
% strokes_topology:
% ------------------------------
% 
%   61×1 struct array with fields:
% 
%     points2D
%     primitive_type
%     primitive_geom
%     mean_pressure
%     length2DFull
%     length3D
%     length2DPrimitive
%     line_group
%     speed
%     accuracy_radius
% ------------------------------
% 
% -------------------------------------------------------------------------
% Output:
% -------------------------------------------------------------------------
% ------------------------------
% strokes_topology:
% ------------------------------
% 
%   num_strokes×1 struct array with fields:
% 
%     points2D
%     primitive_type
%     primitive_geom
%     mean_pressure
%     length2DFull
%     length3D
%     length2DPrimitive
%     line_group
%     speed
%     accuracy_radius
%     poly2d_extended --- added an extension of strokes by accuracy radius
% ------------------------------

function strokes_topology = extendStrokes(strokes_topology, img)
    
    for i = 1:length(strokes_topology)
        % Skip marks:
        strokes_topology(i).poly2d_extended = [];
        
        if strokes_topology(i).primitive_type ~= 0 && ...
           strokes_topology(i).primitive_type ~= 1
            continue;
        end
            
        % Straight:
        if strokes_topology(i).primitive_type == 0
            [begin_s1_x, begin_s1_y, end_s1_x, end_s1_y] = extendIntersectingPolyline(strokes_topology(i));
        end
        
        % Curves:
        if strokes_topology(i).primitive_type == 1
            [begin_s1_x, begin_s1_y, end_s1_x, end_s1_y] = extendCurve(strokes_topology(i), img);
        end
        
        % Assign extended polyline:
        strokes_topology(i).poly2d_extended = [[begin_s1_x; cat(1,strokes_topology(i).points2D.x); end_s1_x],...
                                              [begin_s1_y; cat(1,strokes_topology(i).points2D.y); end_s1_y]];
        strokes_topology(i).poly2d_extended = unique(strokes_topology(i).poly2d_extended, 'rows');
    end

end


function [begin_x, begin_y, end_x, end_y] = extendCurve(sketch_stroke, img)
    coord2D(:,1) = cat(1,sketch_stroke.points2D.x);
    coord2D(:,2)=  cat(1,sketch_stroke.points2D.y);
    %% Begin
    dir = coord2D(1,:)- coord2D(2,:);
    dir = dir/norm(dir);
    begin_x = coord2D(1,1)  + 2*sketch_stroke.accuracy_radius*dir(1);
    begin_y = coord2D(1,2)  + 2*sketch_stroke.accuracy_radius*dir(2);
    
    %% End
    dir = coord2D(end,:)- coord2D(end-1,:);
    dir = dir/norm(dir);
  
    end_x = coord2D(end,1)  + 2*sketch_stroke.accuracy_radius*dir(1);
    end_y = coord2D(end,2)  + 2*sketch_stroke.accuracy_radius*dir(2);      
    
%     figure(15);
%     hold off;
%     imshow(img);
%     hold on;
%     plot(coord2D(:,1),coord2D(:,2),'-*');
%     plot([coord2D(1,1) begin_x],[coord2D(1,2) begin_y],'-o');
%     plot([coord2D(end,1) end_x], [coord2D(end,2) end_y],'o');
end