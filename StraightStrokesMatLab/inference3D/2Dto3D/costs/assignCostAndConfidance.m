function [ind_line_most_probable, max_cost, confidence] = assignCostAndConfidance(cost_full)
%     global USE_FULL
    
%     if ~USE_FULL
        [max_cost, ind_line_most_probable] = max(cost_full);
%     else
%         [max_cost, ind_line_most_probable] = max(cost_full_solution);
%     end
    
    [~, ind_sort] = sort(cost_full, 'descend');
    if length(ind_sort) > 1 
%         confidence = 1 - cost_full(ind_sort(2))/cost_full(ind_sort(1));
        confidence = 1 - cost_full(ind_sort(2))/cost_full(ind_sort(1));
    else
        confidence = 1.0;
    end
           
end