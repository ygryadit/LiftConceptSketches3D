% function [strokes, mask_lines_kept] = computeAverageStrokesSpeed(strokes)
function [strokes] = computeAverageStrokesSpeed(strokes)
%COMPUTEAVERAGESTROKESSPEED 
% 
% Computes the average speed of the stroke.

% strokes_speed = zeros(length(strokes),1);
global datatset_name;

mask_non_defined = false(length(strokes),1);
% inds_lines = find(cat(1,strokes(:).primitive_type) == 0);
% mask_non_defined_lines = false(length(strokes),1);

for  i = 1:length(strokes)
    
    if strcmp(datatset_name, 'OpenSketch')

        % Distances:
        points = [cat(1,strokes(i).points(:).x) cat(1,strokes(i).points(:).y)];
   
        points_start = points(1:(end-1), :);
        points_end   = points(2:end, :);

        distances = sqrt(sum((points_end - points_start).^2, 2));
        distance = sum(distances);

        % Durations:
        times = cat(1,strokes(i).points(:).t);
        times_start = times(1:(end-1), :);
        times_end   = times(2:end, :);

        durations = abs(times_end - times_start);
        duration = sum(durations);


        strokes(i).speed = distance./duration;
    end

    if strcmp(datatset_name, 'AnalyticDrawing')
        velocities = [cat(1,strokes(i).points(:).v)];
        
        num_points = size(velocities);
        
        mean_velocity = sum(velocities)/num_points(1);
        
        strokes(i).speed = mean_velocity;
    end

    if isinf(strokes(i).speed)
        
        mask_non_defined(i) = true;
    end
    
    if (strokes(i).speed < 1e-6)
        
        mask_non_defined(i) = true;
    end
    
    if isnan(strokes(i).speed)
        
        mask_non_defined(i) = true;
    end
    
    if mask_non_defined(i) 
%         mask_non_defined_lines(i) = true;
        strokes(i).primitive_type = -3;
    end
end

% fprintf('computeAverageStrokesSpeed: num_strokes before = %d\n', length(strokes));


% strokes = strokes(~mask_non_defined);
% mask_lines_kept = ~mask_non_defined_lines(inds_lines);


end

