% 25.07.2019 refactor so that 
% 
% strokes has added fileds:
%   Input:
%   strokes:
%       points:
% 
%   Output:
%   strokes_topology:
%       points2D: 
%           copy strokes.points
%       primitive_type: 
%           0 -- line
%           1 -- bezier
%           2 -- mark -- very short line
%       primitive_geom:
%           in case 'is_line' is true frour control poitns [x1 x2 y1 y2]
%           otherwise will be filled in later.
%       mean_pressure:
%           mean pressure of the stroke     
% 
% function [lines, curves, all, lines_mean_pressure] = findAllStrokeLines(strokes, img)
function [strokes_topology] = findAllStrokeLines(strokes)
    
    num_strokes = round(length(strokes)); 
    
    

%     lines [x1 x2 y1 y2 str_num]
    % Some strokes in the end can be removed:
%     fraction_strokes = 1.0;
    
%     lines.coordinates2D = zeros(num_strokes, 4);
%     lines_mean_pressure = zeros(num_strokes, 1);
%     num_lines = 0;
%     num_curves = 0;
%     curve_stroke_nums = [];

    strokes_topology = repmat(struct('points2D',[], ...
                                     'points2DOriginal', [],...
                                     'primitive_type',4, ...
                                     'primitive_geom', [], ...
                                     'mean_pressure', NaN,...
                                     'speed', NaN,...
                                     'line_group', 6,...
                                     'ind_orth_ax', NaN), num_strokes, 1);
    %global filepath_sketch_img;             
    %img = readSketchImg(filepath_sketch_img);
    
    %    imshow(img); hold on;
        
    for i = 1:num_strokes
        stroke = strokes(i).points;
        num_points = length(stroke);
        
        % Assign:
        strokes_topology(i).points2D = stroke;
      
        % Intialise:
        strokes_topology(i).mean_pressure   = strokes(i).mean_pressure;%mean(cat(1,stroke(:).p));
        strokes_topology(i).speed   = strokes(i).speed;%mean(cat(1,stroke(:).p));
        ls = lengthStroke([cat(1, strokes_topology(i).points2D(:).x) cat(1, strokes_topology(i).points2D(:).y)]);
        strokes_topology(i).length2DFull = ls;
        
        
        strokes_topology(i).length3D = [];
        
        strokes_topology(i).length2DPrimitive = NaN;
        
        % Non-defined speed:
        if ~isempty(strokes(i).primitive_type)
            strokes_topology(i).primitive_type = strokes(i).primitive_type;
            continue;
        end
        
        % Mark:
        if num_points <= 2
%             stroke = trimStraightStrokeCoordinatesDrawingField(stroke);
            strokes_topology(i).points2D = stroke;
            
            strokes_topology(i).primitive_geom = [cat(1,stroke(:).x)' cat(1,stroke(:).y)'];
            num_points = length(stroke);
            if num_points == 2
                if ~isempty(strokes_topology(i).primitive_geom)
                    strokes_topology(i).length2DPrimitive = ...
                        lengthStroke([  strokes_topology(i).primitive_geom(1:2)'...
                                        strokes_topology(i).primitive_geom(3:4)']);
                end
                strokes_topology(i).primitive_type = 0;
             
            else
                strokes_topology(i).length2DPrimitive = 0;
                strokes_topology(i).primitive_type = -3;
            end
            continue;
        end
        
        
        
        
        [x,y] = strokeStruct2Vectors(stroke);

        
        %plot(x,y);
% %         

        
        [is_line, line] = isStrokeALine(x,y);
        
        
        if (is_line)
            % Line
%             stroke = trimStraightStrokeCoordinatesDrawingField(stroke);
            strokes_topology(i).points2D = stroke;
            
            strokes_topology(i).primitive_type = 0;

            [x_, y_,~]= trimLineCoordinatesDrawingField(line(1:2), line(3:4));
            strokes_topology(i).primitive_geom = [x_,y_]; 

            %plot(line(:,[1,2]),line(:,[3,4]), ':', 'LineWidth', 2);
            
%            num_lines = num_lines+1;               
%            lines.coordinates2D(num_lines, :) = line(1, 1:4);
%            
%            lines.stroke_nums(num_lines) = i;
            if ~isempty(strokes_topology(i).primitive_geom)
                strokes_topology(i).length2DPrimitive = lengthStroke([strokes_topology(i).primitive_geom(1:2)' strokes_topology(i).primitive_geom(3:4)']);
            end
        else
            % Bezier
            strokes_topology(i).primitive_type = 1;
            
%            num_curves = num_curves + 1;
%            curves{num_curves} = strokes(i).points;
%            curve_stroke_nums = [curve_stroke_nums i];
        end
    
        
    end
% 
%     lines.coordinates2D = lines.coordinates2D(1:num_lines, :);
%     
%     curves_in = curves;
%     clear 'curves';
%     for i = 1:length(curves_in)
%        curves(i).coordinates2D = [cat(1,curves_in{i}.x) cat(1,curves_in{i}.y)];
%        curves(i).length = lengthStroke(curves(i).coordinates2D);
%        curves(i).original_stroke_num = curve_stroke_nums(i);
%     end
end

function [x,y] = strokeStruct2Vectors(stroke)
       x = cat(1,stroke(:).x);
       y = cat(1,stroke(:).y);
end