function [inds_non_assigned_strks,...
        max_costs,...
        confidence_vals,...
        confidence_threshold] = ...
            getNonAssigenedStrokesData(strokes_topology, threshold, inds_non_assigned_strks)
% Find the strokes that are not assigned:


mask = ~[strokes_topology(inds_non_assigned_strks).depth_assigned];
inds_non_assigned_strks = inds_non_assigned_strks(mask);

if ~exist('inds_non_assigned_strks', 'var')
    inds_non_assigned_strks = find(cat(1,  strokes_topology(:).num_candidate_lines) > 0);
end
max_costs = cat(1, strokes_topology(inds_non_assigned_strks).score);
confidence_vals = cat(1, strokes_topology(inds_non_assigned_strks).confidence);

% Keep only strokes with high score values:
max_costs_mask = max_costs > threshold;
max_costs = max_costs(max_costs_mask);
confidence_vals = confidence_vals(max_costs_mask);
inds_non_assigned_strks = inds_non_assigned_strks(max_costs_mask);

if isempty(max_costs)
   confidence_threshold = Inf;
   return;
end

% Sort the strokes according to cost and then according to confidence:  
% [vals, inds] = sortrows([confidence_vals max_costs], 'descend'); 
% confidence_vals = confidence_vals(inds);
% inds_non_assigned_strks = inds_non_assigned_strks(inds);
% max_costs = max_costs(inds);

[inds_non_assigned_strks,...
          max_costs, ...
          confidence_vals] = ...
                sortStrokesScoreConfidence(inds_non_assigned_strks,...
                                    max_costs, ...
                                    confidence_vals);

confidence_threshold = confidence_vals(1);
end