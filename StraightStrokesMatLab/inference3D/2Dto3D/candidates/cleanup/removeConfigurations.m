function [strokes_topology, intersections,strks_cleaned] = removeConfigurations(ind_intrsctng_strk,...
                              ind_candidate_line,...
                              inds_configurations_remove,...
                              strokes_topology,...
                              intersections,...
                              inds_strks_visited)
% 
% if  (ind_intrsctng_strk == 75) & ind_candidate_line == 6
%    disp('check'); 
% end
                          
strks_cleaned =  [inds_strks_visited, ind_intrsctng_strk];                       
% Remove the intersection from confiugrations to remove:
global DISPLAY_INFO
% if DISPLAY_INFO
%    
%     message = 'inds_strks_visited';
%     for i = 1:length(inds_strks_visited)
%         message = sprintf('%s %d', message, inds_strks_visited(i));
%     end
%     displayLog(message);
%     
%     fprintf('removing strk_%d->candaite line: %d, configurations: ',ind_intrsctng_strk,ind_candidate_line);
%     disp(inds_configurations_remove);
%     
% end
try
    candidate_line = strokes_topology(ind_intrsctng_strk).candidate_lines(ind_candidate_line);
catch e
    rethrow(e)
end

candidate_line = remove_references_to_assigned_strokes(...
                    inds_configurations_remove,...
                    candidate_line,...
                    inds_strks_visited);
 

strokes_topology(ind_intrsctng_strk).candidate_lines(ind_candidate_line) = candidate_line;
                          
% Clean up the intersections structure:
% fprintf('ready to clean %d cnd line %d\n', ind_intrsctng_strk, ind_candidate_line);
try
    [intersections, strokes_topology, strks_cleaned] = ...
        cleanCandidateIntersections(inds_configurations_remove,...
                                    strokes_topology,...
                                    ind_intrsctng_strk,...
                                    ind_candidate_line,...
                                    intersections,...
                                    inds_strks_visited);
catch e
    rethrow(e);
end
% fprintf('done cleaning %d cnd line %d\n', ind_intrsctng_strk, ind_candidate_line);
    
% Should remove configurations and reindex all the indices of configurationd
% in intersectins in remaining configurations:
inds_configurations_all = 1:length(strokes_topology(ind_intrsctng_strk).candidate_lines(ind_candidate_line).configurations);
inds_keep = setdiff(inds_configurations_all, inds_configurations_remove);

strokes_topology(ind_intrsctng_strk).candidate_lines(ind_candidate_line).configurations = ...
    strokes_topology(ind_intrsctng_strk).candidate_lines(ind_candidate_line).configurations(inds_keep);

% if isempty(inds_keep)
%     [strokes_topology, intersections] = ...
%         removeCandidateLine(strokes_topology, ...
%                             ind_intrsctng_strk,...
%                             ind_candidate_line ,...
%                             intersections);
% else
try
    intersections = reindexConfigurations(inds_keep, strokes_topology, ind_intrsctng_strk, ind_candidate_line, intersections );
catch e
    rethrow e
end
% end

%% Update scores;
strokes_topology(ind_intrsctng_strk) = computeMaxScoreStrokeBranches(strokes_topology(ind_intrsctng_strk));
%% Update list joint strokes;
strokes_topology(ind_intrsctng_strk).candidate_lines = ...
    updateListJointStrokes(strokes_topology(ind_intrsctng_strk).candidate_lines);

end


function candidate_line = remove_references_to_assigned_strokes(inds_configurations_remove,...
                candidate_line,...
                inds_strks_visited)

    for i = 1:length(inds_configurations_remove)
        if isempty( candidate_line.configurations)
            continue;
        end
        cnfgrtn = candidate_line.configurations(inds_configurations_remove(i));

        mask_intrsctns_kp = ~ismember(cnfgrtn.inds_jnts_strks, inds_strks_visited);
        cnfgrtn.inds_intrsctns__mult_cnddts = ...
            cnfgrtn.inds_intrsctns__mult_cnddts(mask_intrsctns_kp);
        cnfgrtn.inds_intrsctns__mult_cnddts_ind = ...
            cnfgrtn.inds_intrsctns__mult_cnddts_ind(mask_intrsctns_kp);
        cnfgrtn.inds_jnts_strks = ...
            cnfgrtn.inds_jnts_strks(mask_intrsctns_kp);


    %     mask_intrsctns_kp = ~ismember(cnfgrtn.inds_intrsctns, inds_strks_visited);
        inds_intrsctns_keep = [cnfgrtn.inds_intrsctns__assigned cnfgrtn.inds_intrsctns__mult_cnddts];
        ind_keep = ismember( cnfgrtn.inds_intrsctns, inds_intrsctns_keep);
        cnfgrtn.inds_intrsctns = inds_intrsctns_keep;
    %         cnfgrtn.inds_intrsctns(mask_intrsctns_kp);
        cnfgrtn.p_intrsctns_dists = cnfgrtn.p_intrsctns_dists(ind_keep); 
    %         cnfgrtn.p_intrsctns_dists(mask_intrsctns_kp); 

        candidate_line.configurations(inds_configurations_remove(i)) = cnfgrtn;
    end   

end