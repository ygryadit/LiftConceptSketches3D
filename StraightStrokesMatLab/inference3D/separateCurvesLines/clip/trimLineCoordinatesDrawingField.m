% function strokes_topology = trimStraightStrokesCoordinatesDrawingField(strokes_topology)
% 
% Description:
% In OpenSketch strokes were recorded outside the boundary of the sketching
% area if the stroke started inside the sketchign area. This funstion trims
% straight strokes to lie within the original sketching area.

function [x_vals, y_vals, inside] = ...
                trimLineCoordinatesDrawingField( ...
                    x_vals,...
                    y_vals,...
                    debug)
    global sketch_height;
    
    % Region of clipping:
    min_b = 0;
    max_b = sketch_height;

    
    if ~exist('debug', 'var')
        debug = false;
    end

    %% Determine strokes directions:
    
    grows_x = (x_vals(end) - x_vals(1)) > 0;
    grows_y = (y_vals(end) - y_vals(1)) > 0;
    
    %% x 0
    inds_out = find(x_vals < min_b);
    
    [x_vals, y_vals] = ...
            updateOneClipCoordinateMin(x_vals,...
                                    y_vals,...
                                    inds_out,...
                                    grows_x, ...
                                    'x',...
                                    min_b);

    
    %% x - height
    inds_out = find(x_vals > max_b);
    
%     if ~isempty(inds_out)
%         if grows_x
%             vert_num = inds_out(1);
%             seg = [x_vals(vert_num) y_vals(vert_num); ...
%                    x_vals(vert_num-1) y_vals(vert_num-1)];  
%         else
%             vert_num = inds_out(end);
%             seg = [x_vals(vert_num)  y_vals(vert_num); ...
%                    x_vals(vert_num+1) y_vals(vert_num+1)];
%         end
%         
%         y_vals(vert_num) = interpolate_x(seg, max_b);
%         x_vals(vert_num) = max_b;
%     end
        [x_vals, y_vals] = ...
            updateOneClipCoordinateMax(x_vals,...
                                    y_vals,...
                                    inds_out,...
                                    grows_x, ...
                                    'x',...
                                    max_b);
    
    %% y 0 
    inds_out = find(y_vals < min_b);
       [x_vals, y_vals] = ...
            updateOneClipCoordinateMin(x_vals,...
                                    y_vals,...
                                    inds_out,...
                                    grows_y, ...
                                    'y',...
                                    min_b);
    
    
    %% y - height
    inds_out = find(y_vals > max_b);
    
    [x_vals, y_vals] = ...
            updateOneClipCoordinateMax(x_vals,...
                                    y_vals,...
                                    inds_out,...
                                    grows_y, ...
                                    'y',...
                                    max_b);
    
    %% Define points inside:
    
    inside = find(  x_vals >= min_b & ...
                    x_vals <= max_b & ...
                    y_vals >= min_b & ...
                    y_vals <= max_b);

    % Keep only strokes inside the sketching area:
%     stroke = stroke(inside);
    
    x_vals = x_vals(inside);
    y_vals = y_vals(inside);

    if debug
        plot(x_vals, y_vals);
    end


end






