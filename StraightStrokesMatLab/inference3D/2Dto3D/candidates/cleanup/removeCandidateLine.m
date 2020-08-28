function [strokes_topology, intersections] = ...
        removeCandidateLine(strokes_topology, ...
                            cur_stroke,...
                            inds_cnddt_lns_rmv ,...
                            intersections)
ind_stroke = cur_stroke.ind;                        
candidate_lines = strokes_topology(ind_stroke).candidate_lines;

% if ind_stroke == 33 && ismember(5, inds_cnddt_lns_rmv)
%    disp('check');
% end

% Indices of candiate lines to keep:
num_candidate_lines = length(candidate_lines);
inds_keep = setdiff(1:num_candidate_lines, inds_cnddt_lns_rmv);

%% Remove configurations from the candidate line
for i = inds_cnddt_lns_rmv
    cnfgrtns = strokes_topology(cur_stroke.ind).candidate_lines(i).configurations;
    if ~isempty(cnfgrtns)
        try
        [strokes_topology,...
          intersections] = removeConfigurations(cur_stroke.ind,...
                                          i,...
                                          1:length(cnfgrtns),...
                                          strokes_topology,...
                                          intersections,...
                                          []);                    
        catch e
            rethrow(e);
        end
    end
end
                

%% Remove references to candiates lines to remove 
% Go over all the intersections with multiple candidates


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

for ind_intrsctn = reshape(cur_stroke.inds_intrsctns_eval_mltpl_cnddts, 1,[])
    ind_pair = find(intersections(ind_intrsctn).strokes_indices == ind_stroke);
    
    if ~isfield(intersections(ind_intrsctn), 'cnddts3D')
            %The sroke did not have any candidate lines.
            continue;
    end
    
    for intrsctn_vrsn = 1:length(intersections(ind_intrsctn).cnddts3D)
        intrsctn = intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn);     
        intrsctn.cnddt_lns{ind_pair} = setdiff(intrsctn.cnddt_lns{ind_pair}, inds_cnddt_lns_rmv);
        intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn) = intrsctn;
    end
end


% % Remove reference to candidate line
% for i = inds_cnddt_lns_rmv
%     % Go over the configurations
%     for j = 1:length(candidate_lines(i).configurations)
%         configuration = candidate_lines(i).configurations(j);
%         % Find itersections with multiple candidates:
%         inds_intrsctns__mult_cnddts = ...
%             configuration.inds_intrsctns__mult_cnddts;
%         inds_intrsctns__mult_cnddts_ind = ...
%             configuration.inds_intrsctns__mult_cnddts_ind;
%         
%         % Perform reindexing
%         for ijk  = 1:length(inds_intrsctns__mult_cnddts)
%             ind_intrsctn = inds_intrsctns__mult_cnddts(ijk);
%             intrsctn_vrsn = inds_intrsctns__mult_cnddts_ind(ijk);
%             intrsctn = intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn);
%             
%             ind_pair = find(intersections(ind_intrsctn).strokes_indices == ind_stroke);
%             
%             intrsctn.cnddt_lns{ind_pair} = setdiff(intrsctn.cnddt_lns{ind_pair}, i);
%             
%             
%             intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn) = intrsctn;
%         end
%     end
% end



%% Remove the candidate line:
candidate_lines = candidate_lines(inds_keep);

%% Reindex the remaining candidate lines:
% 	for each intersection with multiple candidates change the number of the candidate line
num_candidate_lines = length(candidate_lines);
inds_new = 1:num_candidate_lines;

for ind_intrsctn = reshape(cur_stroke.inds_intrsctns_eval_mltpl_cnddts, 1,[])
    ind_pair = find(intersections(ind_intrsctn).strokes_indices == ind_stroke);
    
    if ~isfield(intersections(ind_intrsctn), 'cnddts3D')
       %The sroke did not have any candidate lines.
       continue;
    end
    
    for intrsctn_vrsn = 1:length(intersections(ind_intrsctn).cnddts3D)
        intrsctn = intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn);     
        
        mask_update = ismember(intrsctn.cnddt_lns{ind_pair}, inds_keep);
        mask_new_vals = ismember(inds_keep, intrsctn.cnddt_lns{ind_pair});
        
        intrsctn.cnddt_lns{ind_pair}(mask_update) = inds_new(mask_new_vals);
        
        intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn) = intrsctn;
    end
end




% for i = 1:num_candidate_lines
%     % Go over the configurations
%     for j = 1:length(candidate_lines(i).configurations)
%         configuration = candidate_lines(i).configurations(j);
%         % Find itersections with multiple candidates:
%         inds_intrsctns__mult_cnddts = ...
%             configuration.inds_intrsctns__mult_cnddts;
%         inds_intrsctns__mult_cnddts_ind = ...
%             configuration.inds_intrsctns__mult_cnddts_ind;
%         
%         % Perform reindexing
%         for ijk  = 1:length(inds_intrsctns__mult_cnddts)
%             ind_intrsctn = inds_intrsctns__mult_cnddts(ijk);
%             intrsctn_vrsn = inds_intrsctns__mult_cnddts_ind(ijk);
%             intrsctn = intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn);
%             
%             ind_pair = find(intersections(ind_intrsctn).strokes_indices == ind_stroke);
%             
%             intrsctn.cnddt_lns{ind_pair}(intrsctn.cnddt_lns{ind_pair} == inds_keep(i)) = i;
%             
%             intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn) = intrsctn;
%         end
%     end
% end

if isempty(candidate_lines)
    candidate_lines = [];
end

strokes_topology(ind_stroke).candidate_lines = candidate_lines;

strokes_topology(ind_stroke).num_candidate_lines = length(candidate_lines);
                                    
                                    
end