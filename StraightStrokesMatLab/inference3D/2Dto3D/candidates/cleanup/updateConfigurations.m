
function [strokes_topology,...
          intersections] = updateConfigurations(ind_stroke,...
                                                ind_candidate_line,...
                                                inds_configurations_update,...
                                                ind_intrsctn_update,...
                                                strokes_topology,...
                                                intersections)
        
%      strokes_topology = removeConfigurations(ind_stroke,...
%                                   ind_candidate_line,...
%                                   inds_configurations_update,...
%                                   strokes_topology);
                              
%      configurations = strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations;                        
%      num_configs = length(configurations);
%      for i = 1:num_configs
%         configurations(i).inds_intrsctns__assigned(end+1) = ind_intrsctn_update;    
%      end
% if  ind_intrsctn_update == 416
%    disp('check'); 
% end   
%                 
% Update    
for i = 1:length(inds_configurations_update)
    ind_cnfgrtn_update = inds_configurations_update(i);
    if isempty(strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations)
        warning('Safer to do a cleanup, but difficult so trying to get away with this')
        continue;
    end
    try
    configuration = strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations(ind_cnfgrtn_update);
    catch e
        fprintf('ind_stroke = %d, ind_candidate_line = %d, ind_cnfgrtn_update = %d\n ', ...
                 ind_stroke, ind_candidate_line, ind_cnfgrtn_update)
        rethrow(e);     
    end
        
    mask_intrsctn = ismember(configuration.inds_intrsctns__mult_cnddts, ind_intrsctn_update);
    configuration.inds_intrsctns__mult_cnddts = configuration.inds_intrsctns__mult_cnddts(~mask_intrsctn);
    configuration.inds_intrsctns__mult_cnddts_ind = configuration.inds_intrsctns__mult_cnddts_ind(~mask_intrsctn);
    configuration.inds_jnts_strks = configuration.inds_jnts_strks(~mask_intrsctn);
    configuration.inds_intrsctns__assigned(end+1) = ind_intrsctn_update;
    strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations(ind_cnfgrtn_update) = configuration;   
    
    
end





% Remove the rest   
if intersections(ind_intrsctn_update).is_active == 1
    if ~isempty(strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations)
%         inds_configurations_remove = setdiff(1:length( strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations),...
%                                             inds_configurations_update);
%         if (ind_stroke == 36)
%             disp('check'); 
%         end

        strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations = ...
           strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations(inds_configurations_update);
       intersections = reindexConfigurations(inds_configurations_update, strokes_topology, ind_stroke, ind_candidate_line, intersections );
    end
end


%% Update scores;

strokes_topology(ind_stroke) = ...
    computeMaxScoreStrokeBranches(strokes_topology(ind_stroke));

%% Update list joint strokes;
strokes_topology(ind_stroke).candidate_lines = ...
    updateListJointStrokes(strokes_topology(ind_stroke).candidate_lines);

end