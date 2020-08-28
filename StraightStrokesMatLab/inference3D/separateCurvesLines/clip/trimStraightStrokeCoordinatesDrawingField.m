% function strokes_topology = trimStraightStrokesCoordinatesDrawingField(strokes_topology)
% 
% Description:
% In OpenSketch strokes were recorded outside the boundary of the sketching
% area if the stroke started inside the sketchign area. This funstion trims
% straight strokes to lie within the original sketching area.

function stroke = trimStraightStrokeCoordinatesDrawingField(stroke)
    global sketch_height;
    
    % Region of clipping:
%     min_b = 0;
%     max_b = sketch_height;

    debug = false;
    
    if debug
        h1 = figure(1);
        hold on;
        axis equal;
    end

%     % Evaluate which strokes points are inside the sketching area:
    x_vals = cat(1,stroke.x);
    y_vals = cat(1,stroke.y);
    
    [~, ~, inside] = ...
        trimLineCoordinatesDrawingField(x_vals, y_vals, debug);

    % Keep only strokes inside the sketching area:
    stroke = stroke(inside);



end

function coord = interpolate_x(seg, val_new)
    y2 = seg(2,2); 
    y1 = seg(1,2);
    
    x2 = seg(2,1);
    x1 = seg(1,1);
    
    dy = y2 - y1;
    y_grows = dy > 0;
    dx = x2 - x1;
    
    r = dy/dx;
    
    if y_grows
        coord = r*(val_new - x1);
    else
        coord = r*(val_new - x2);
    end

end


function coord = interpolate_y(seg, val_new)
    y2 = seg(2,2); 
    y1 = seg(1,2);
    
    x2 = seg(2,1);
    x1 = seg(1,1);
    
    dx = x2 - x1;
    dy = y2 - y1;
    x_grows = dx > 0;
    
    
    r = dx/dy;
    
    if x_grows
        coord = r*(val_new - y1);
    else
        coord = r*(val_new - y2);
    end

end

