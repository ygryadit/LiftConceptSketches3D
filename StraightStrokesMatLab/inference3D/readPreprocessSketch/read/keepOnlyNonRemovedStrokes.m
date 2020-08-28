function sketch_out = keepOnlyNonRemovedStrokes(sketch)
% Input:
%   sketch: struct with fields:
% 
%       canvas: struct with fields:
%            height: sketch height
%            pen_width: 1.5000
%            width: sketch width
%         
%       strokes: struct with fields:
%           is_removed
%           points

% Output:
%   sketch: struct with fields:
% 
%       canvas: struct with fields:
%            height: sketch height
%            pen_width: 1.5000
%            width: sketch width
%         
%       strokes: struct with fields:
%           points
  
    global IPad;
    IPad = false;
    if isfield(sketch, 'canvas')
        sketch_out.canvas = sketch.canvas;
        
    else
        sketch_out.canvas.height = NaN;
        sketch_out.canvas.width = NaN;
        sketch_out.canvas.pen_width = 1.5;
    end
    
   if ~isfield(sketch_out.canvas, 'height')
       sketch_out.canvas.height = sketch_out.canvas.width;
       sketch_out.canvas.pen_width = 1.5;
   end
    
    num_strokes_in = length(sketch.strokes);
    num_strokes_out = 0;
    
    
    for i = 1:num_strokes_in
        if (iscell(sketch.strokes))
            stroke = sketch.strokes{i};
        else
            stroke = sketch.strokes(i);
        end
        
        if isfield(stroke, 'is_removed') 
            if ( ~stroke.is_removed && isfield(stroke, 'points'))
                num_strokes_out = num_strokes_out+1;
                sketch_out.strokes(num_strokes_out).points = ...
                    stroke.points;
            end
        else
            IPad = true;
            num_strokes_out = num_strokes_out+1;
            sketch_out.strokes(num_strokes_out).points = ...
                stroke.points2D;
        end
    end 
end