function [distances_begin, distances_end] = computeIntersectionsCurveEndPointsDistance(intrsctns_seg_nums, stroke)
    distances_begin = zeros(length(intrsctns_seg_nums),1);
    distances_end   = zeros(length(intrsctns_seg_nums),1);

    
    
    for i = 1:length(intrsctns_seg_nums)

        distances_begin =  computeStrokeLengthBetweenSegPositions(stroke, intrsctns_seg_nums(1), intrsctns_seg_nums(i));
        distances_end =  computeStrokeLengthBetweenSegPositions(stroke, intrsctns_seg_nums(i), intrsctns_seg_nums(end));
 
    end
end