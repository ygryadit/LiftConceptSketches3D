function inds_configurations = getWinningConfigurations(p_joint_configurations, candidate_line)
       [ind_configuration_most_probable, max_cost, confidence] = ...
            assignCostAndConfidance(p_joint_configurations);
        
        if doNotAssignDepthValueToLine(max_cost,  confidence, max_cost, length(candidate_line.configurations), candidate_line)
            inds_configurations = 1:length(p_joint_configurations);
        else
            inds_configurations = ind_configuration_most_probable;   
        end
end