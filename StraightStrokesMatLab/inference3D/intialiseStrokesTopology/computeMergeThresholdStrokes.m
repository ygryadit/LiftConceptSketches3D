% function [strokes_thresholds] = computeMergeThresholdStrokes(strokes_speed, MIN_MERGE_THR, MAX_MERGE_THR)
function [strokes_thresholds] = computeMergeThresholdStrokes(strokes_speed, scale, MIN_MERGE_THR, MAX_MERGE_THR)
%COMPUTEMERGETHRESHOLDSTROKES 
min_s_sketch = min(strokes_speed);
max_s_skecth = max(strokes_speed);

%     [MIN_MERGE_THR,MAX_MERGE_THR] = computeGroupingThreshold(sketch);
     MIN_MERGE_THR = 3*scale; %2*line_width
     MAX_MERGE_THR = 10*scale;


mean_stroke_spead = mean(strokes_speed);
std_spead = std(strokes_speed);

%% Values from our dataset:
min_s = 0.3929; 
max_s = 1.7940e+03;
% mean_stroke_spead - 3*std_spead;
% max_s = mean_stroke_spead + 3*std_spead;


strokes_thresholds = MIN_MERGE_THR + (MAX_MERGE_THR - MIN_MERGE_THR)./(max_s - min_s).*(strokes_speed);


end

