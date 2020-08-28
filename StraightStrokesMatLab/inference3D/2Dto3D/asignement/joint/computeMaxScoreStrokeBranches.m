% function computeMaxScoreStrokeBranches()
% 
% Description:
%   Goes over candidate lines and configurations and selects max.
function stroke = computeMaxScoreStrokeBranches(stroke)
    
    for i =1:length(stroke.candidate_lines)
        stroke.candidate_lines(i).max_cost = max(cat(1, stroke.candidate_lines(i).configurations(:).p_full_joint));
    end
    vals = cat(1, stroke.candidate_lines(:).max_cost);
   [~,  stroke.score, ...
        stroke.confidence] = ...
                assignCostAndConfidance(vals);
    if isempty(stroke.score)
       stroke.score  = 0;
       stroke.candidate_lines = [];       
    end        
   
end