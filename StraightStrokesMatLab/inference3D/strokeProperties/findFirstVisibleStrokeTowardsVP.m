%   probabilitiesVP:
%       num_straight_strokes x 4 (3 vanishign points and the rest)
%       votes normilised by the length of each stroke

function ind_stroke = ...
    findFirstVisibleStrokeTowardsVP(...
                probabilitiesVP,...
                strokes_topology)

    ind_stroke = 1; % Index of the stroke in the strokes topology
    lj = 1; % Index of the stroke among consecutive straight strokes
        
    [~, lines_group] = max(probabilitiesVP,[],2);  

    while strokes_topology(ind_stroke).primitive_type ~= 0 || ...
          lines_group(lj) == 4 || ...
          (strokes_topology(ind_stroke).mean_pressure < 0.05) 

        if (lines_group(lj) == 4)
           lj = lj + 1; 
        end
         ind_stroke = ind_stroke+1;
    end

end