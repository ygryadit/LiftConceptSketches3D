
function val =  doNotAssignDepthValueToLine(max_cost, confidence, max_cost_stroke, num_lines, candidate_lines)
global confidence_threshold;
global max_cost_threshold;
global max_cost_absolute;

global ASSIGN_BEST;
if ASSIGN_BEST
    val = false;
    return;
end

global ASSIGN_HIGH_SCORE;
if num_lines==1
     disp('');
end

% global UPDATE_STROKES;
% if UPDATE_STROKES
%     val = ((max_cost < 0.01) || ...
%          (max_cost < 0.5 && confidence < 0.2));        
% else
%     val = ((max_cost < 0.01) ||...
%         ( confidence < 0.2));
% end

% assign_depth = (max_cost > 0.75) |  (max_cost > 0.5 & confidence > 0.2);
% assign_depth = (max_cost > 0.75) |  (max_cost > 0.6 & confidence > confidence_threshold);
% assign_depth = (max_cost > 0.85) |  (max_cost > max_cost_threshold & confidence > confidence_threshold);

% assign_depth = ((max_cost > 0.95) | (max_cost > max_cost_threshold) & (confidence > confidence_threshold));
% assign_depth = ((max_cost > 0.85) | (max_cost > max_cost_threshold) & (confidence > confidence_threshold));
if num_lines==1
    assign_depth = true;
else
%     assign_depth = ( (max_cost > 0.98) | ((max_cost > max_cost_absolute) & ( confidence > confidence_threshold2) ) | (max_cost > max_cost_threshold) & (confidence > confidence_threshold));
%       conf_thr = confidence_threshold/0.5*(1 - max_cost);
%       conf_thr = confidence_threshold*sqrt(1.0 - max_cost);
      
%       assign_depth =  (confidence > conf_thr) & (max_cost > max_cost_threshold);
%     assign_depth = ( (max_cost >= 0.98) | (max_cost >= max_cost_threshold) & (confidence >= confidence_threshold));
% assign_depth = ( (max_cost >= max_cost_threshold) & (confidence >= confidence_threshold));
    if ASSIGN_HIGH_SCORE
        assign_depth = ( (max_cost >= max_cost_absolute) | ...
                         (  (max_cost >= max_cost_threshold) &...
                            (confidence >= confidence_threshold) & ...
                            (max_cost_stroke >= max_cost_threshold) ) );
    else
        assign_depth = (  (max_cost >= max_cost_threshold) &...
                            (confidence >= confidence_threshold) & ...
                            (max_cost_stroke >= max_cost_threshold) );
    end
end


% if confidence < confidence_threshold
%    %Check that the solution within confidence radius are similar:
%    %The angle is less than 30 degrees
%    %The distance is less than 1/3 of the length of strokes with minimum 3d
%    %length.
%    thr_cost = max_cost - 0.5*confidence_threshold*max_cost;
%    
%    mask_best = (p_joint_cnddt_lns == max(p_joint_cnddt_lns));
%    
%    mask_simlar_cost = (p_joint_cnddt_lns > thr_cost) & ~mask_best;
%    
%    best_line = candidate_lines(mask_best);
%    candidate_lines = candidate_lines(mask_simlar_cost);
%    
%    
%    % Compare:
%    is_angular_similar = true;
%    is_spatially_near = true;
%    for i = 1:length(candidate_lines)
%       % Angular:      
%       is_angular_similar = is_angular_similar & abs(dot(candidate_lines(i).dir, best_line.dir)) < cos(pi/6);
%       if ~is_angular_similar 
%           % Not similar, do not assign depth:
%           val = true;
%           return;
%       end
%       
%       %Spatial
%       [dist,~] = distLinSeg(candidate_lines(i).coordinates3D_prior(1:3),...
%                                     candidate_lines(i).coordinates3D_prior(4:6),...
%                                     best_line.coordinates3D_prior(1:3),...
%                                     best_line.coordinates3D_prior(4:6));
%       is_spatially_near = is_spatially_near && (dist < 0.3*min(best_line.length3D,  candidate_lines(i).length3D));
%       
%       if ~is_spatially_near 
%           % Not similar, do not assign depth:
%           val = true;
%           return;
%       end
%       
%    end
% end

% assign_depth = ((max_cost > max_cost_threshold) & (confidence > confidence_threshold));
val = ~(assign_depth);