function intersections = reindexConfigurations(inds_keep, strokes_topology, ind_stroke, ind_candidate_line, intersections )
    
    % I. Find all the intersections with multiple candisates and their
    % versions.
    cur_stroke = strokes_topology(ind_stroke);
    cur_stroke.ind = ind_stroke;
   
      UP_TO_LAST = true;
    [cur_stroke.inds_intrsctns_eval,...
     cur_stroke.inds_intrsctns_eval_actv,...
     cur_stroke.inds_intrsctns_eval_mltpl_cnddts,...
     cur_stroke.inds_intrsctng_strks_eval,...
     cur_stroke.inds_intrsctng_strks_eval_actv,...
     cur_stroke.inds_intrsctng_strks_eval_mltpl_cnddts] = ...
        returnIndicesNodesTypes(cur_stroke, ...
                            cat(1, strokes_topology(:).depth_assigned),...
                                        intersections,...
                                        UP_TO_LAST);
    
     inds_remnng_cnfgrts = 1:length(inds_keep);
    for j = 1:length(cur_stroke.inds_intrsctns_eval_mltpl_cnddts)
        ind_inter = cur_stroke.inds_intrsctns_eval_mltpl_cnddts(j);
        
        if ~isfield(intersections(ind_inter), 'cnddts3D')
            %The sroke did not have any candidate lines.
            continue;
        end
        
        for ind_cnddts_ind = 1:length(intersections(ind_inter).cnddts3D)
            ind_pair = find(intersections(ind_inter).strokes_indices == ind_stroke);
        
            inds_cnfgrtns_cell_array = intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnfgrtns{ind_pair};
            ijk = find(intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnddt_lns{ind_pair} == ind_candidate_line);
            if isempty(ijk)
                continue;
            end
            inds_cnfgrtns_before = inds_cnfgrtns_cell_array{ijk};
            mask = ismember(inds_keep, inds_cnfgrtns_before);
            inds_cnfgrtns_after = inds_remnng_cnfgrts(mask);            
            intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnfgrtns{ind_pair}{ijk} = inds_cnfgrtns_after;        
        end
        
    end
    
%     inds_remnng_cnfgrts = 1:length(inds_keep);
%     IM = [];
%     j = 0;
%     candidate_line = strokes_topology(ind_stroke).candidate_lines(ind_candidate_line);
%     for i = inds_remnng_cnfgrts %over all the configurations
%         configuration = candidate_line.configurations(i);
%         for ii = 1:length(configuration.inds_intrsctns__mult_cnddts)
%             j = j+1;
%             IM(j, 1) = configuration.inds_intrsctns__mult_cnddts(ii);
%             IM(j, 2) = configuration.inds_intrsctns__mult_cnddts_ind(ii);
%         end
%     end 
%     
%     IM = unique(IM, 'rows');
    
%     for j = 1:size(IM,1)
%         ind_inter      = IM(j,1);
%         ind_cnddts_ind = IM(j,2);
%         
%         ind_pair = find(intersections(ind_inter).strokes_indices == ind_stroke);
%         
%         inds_cnfgrtns_cell_array = intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnfgrtns{ind_pair};
%         ijk = find(intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnddt_lns{ind_pair} == ind_candidate_line);
%         
%         inds_cnfgrtns_before = inds_cnfgrtns_cell_array{ijk};
%         mask = ismember(inds_keep, inds_cnfgrtns_before);
%         inds_cnfgrtns_after = inds_remnng_cnfgrts(mask);            
%         intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnfgrtns{ind_pair}{ijk} = inds_cnfgrtns_after;
%         
%     end
    
    
    
    
%     inds_remnng_cnfgrts = 1:length(inds_keep);
%     for i = inds_remnng_cnfgrts %over all the configurations
% 
%         configuration = strokes_topology(ind_stroke).candidate_lines(ind_candidate_line).configurations(i);
% 
%         for j = 1:length(configuration.inds_intrsctns__mult_cnddts)
%             ind_inter      = configuration.inds_intrsctns__mult_cnddts(j);
%             ind_cnddts_ind = configuration.inds_intrsctns__mult_cnddts_ind(j);
% 
%             ind_pair = find(intersections(ind_inter).strokes_indices == ind_stroke);
% 
%             inds_cnfgrtns_cell_array = intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnfgrtns{ind_pair};
%             ijk = find(intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnddt_lns{ind_pair} == ind_candidate_line);
%             
% %             for ijk = 1:length(inds_cnfgrtns_cell_array)
%             inds_cnfgrtns = inds_cnfgrtns_cell_array{ijk}
%             mask = ismember(inds_keep, inds_cnfgrtns);
%             inds_cnfgrtns = inds_remnng_cnfgrts(mask)            
%             intersections(ind_inter).cnddts3D(ind_cnddts_ind).cnfgrtns{ind_pair}{ijk} = inds_cnfgrtns;
% %             end
%         end
%     end

end
