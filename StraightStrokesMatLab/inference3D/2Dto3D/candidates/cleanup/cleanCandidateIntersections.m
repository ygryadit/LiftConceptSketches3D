% function [intersections,...
%           strokes_topology] = ...
%               cleanCandidateIntersections(inds_configurations_remove,...
%                                           strokes_topology,...
%                                           ind_stroke,...
%                                           ind_candidate_line,...
%                                           intersections)
% 
% Description:
%   Function is called when the configuration is removed.
%   Find the intersection that had this candidate lines.
%   It either removes configuration or removes the whole candidate
%   intersection and updates the associated structures.

function [intersections, strokes_topology, strks_cleaned] = cleanCandidateIntersections(inds_configurations_remove,...
                                     strokes_topology,...
                                     ind_stroke,...
                                     ind_candidate_line,...
                                     intersections,...
                                     ind_stroke_resolved)

strks_cleaned = [ind_stroke_resolved, ind_stroke];                             
% I. Find all the intersections with multiple candidates and their
% versions that should be removed.
IM = [];
j = 0;
candidate_line = strokes_topology(ind_stroke).candidate_lines(ind_candidate_line);
for i = inds_configurations_remove %over all the configurations to remove
    if isempty( candidate_line.configurations)
        continue;
    end
    configuration = candidate_line.configurations(i);
    for ii = 1:length(configuration.inds_intrsctns__mult_cnddts)
        j = j+1;
        IM(j, 1) = configuration.inds_intrsctns__mult_cnddts(ii);
        IM(j, 2) = configuration.inds_intrsctns__mult_cnddts_ind(ii);
    end
end 

IM = unique(IM, 'rows');
    
for j = 1:size(IM,1)

    ind_intrsctn    = IM(j,1);
    intrsctn_vrsn   = IM(j,2);
    %fprintf('ind_intrsctn_clean = %d version %d\n', ind_intrsctn, intrsctn_vrsn);
    
%     if ind_intrsctn == 80 & intrsctn_vrsn == 4
%         disp('check');
%     end
    
    intrsctn = intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn);

    mask_pair = intersections(ind_intrsctn).strokes_indices == ind_stroke;

    
    [intrsctn, inds_cnfgrtns] = cleanReferencesToRemovedConfigurations( intrsctn, ...%intersection structure version to edit
                                                mask_pair, ...%index of the intersections in a pair of intersecting intersections
                                                ind_candidate_line,... 
                                                inds_configurations_remove);
                                            
    
    intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn) = intrsctn;
    
    
    if isempty(inds_cnfgrtns) && (length(intrsctn.cnddt_lns{mask_pair}) == 0)
        
        % If there are no configuration left => the intersections
        % should be removed.

        % Remove associated configurations lines:
        ind_intrsctng_strk = intersections(ind_intrsctn).strokes_indices(~mask_pair);

        for iii = 1:length(intrsctn.cnddt_lns{~mask_pair})

            ind_candidate_line_ = intrsctn.cnddt_lns{~mask_pair}(iii);
            inds_configurations_remove_ = intrsctn.cnfgrtns{~mask_pair}{iii};
            
            try
            [strokes_topology, intersections, strks_cleaned_] = ...
                removeConfigurations( ind_intrsctng_strk,...
                                      ind_candidate_line_,...
                                      inds_configurations_remove_,...
                                      strokes_topology,...
                                      intersections,...
                                      [ind_stroke_resolved]);%[ind_stroke_resolved ind_stroke]);
            
            catch e
                rethrow(e)
            end
            strks_cleaned = unique([strks_cleaned, strks_cleaned_]);
        end

        
        intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn).cnfgrtns{~mask_pair} = {};
        intersections(ind_intrsctn).cnddts3D(intrsctn_vrsn).cnddt_lns{~mask_pair} =[];

    end
%     strokes_topology(ind_stroke).candidate_lines(ind_candidate_line) = ...
%         removeIntersectionFromCandidateLines(ind_intrsctn, intrsctn_vrsn, candidate_line); 
end
end                             

function candidate_line = removeIntersectionFromCandidateLines(ind_intrsctn, intrsctn_vrsn, candidate_line)
    num_configuration = length(candidate_line.configurations);
    
    for i = 1:num_configuration
        configuration = candidate_line.configurations(i);
        
        mask_delete = ( (configuration.inds_intrsctns__mult_cnddts == ind_intrsctn) & ...
                        (configuration.inds_intrsctns__mult_cnddts_ind == intrsctn_vrsn) ) ;
        
        if ~sum(mask_delete)
            continue;
        end
                    
        mask_keep = ~mask_delete;
        
                        
        configuration.inds_intrsctns__mult_cnddts = ...
            configuration.inds_intrsctns__mult_cnddts(mask_keep);
        
        configuration.inds_jnts_strks = ...
            configuration.inds_jnts_strks(mask_keep);
        
        
        
        mask_keep = configuration.inds_intrsctns ~= ind_intrsctn;
        
        configuration.p_intrsctns_dists = ...
            configuration.p_intrsctns_dists(mask_keep); 
        
        
        
        configuration.inds_intrsctns__mult_cnddts = ...
            setdiff(configuration.inds_intrsctns__mult_cnddts, ind_intrsctn);
        
         configuration.inds_intrsctns = ...
            setdiff(configuration.inds_intrsctns, ind_intrsctn);
        
        candidate_line.configurations(i) = configuration;
    end
     
end



function [intrsctn,... %cleanedup intersection
            inds_cnfgrtns] = ... % remaingin configurations
                cleanReferencesToRemovedConfigurations( intrsctn, ...%intersection structure version to edit
                                                mask_pair, ...%index of the intersections in a pair of intersecting intersections
                                                ind_candidate_line,... 
                                                inds_configurations_remove) %indices configurations to remove
    
    ind_cnddt_ln = find(intrsctn.cnddt_lns{mask_pair} == ind_candidate_line);
    
    if isempty(ind_cnddt_ln)
        inds_cnfgrtns = [];
        return;
    end
                                            
    
    
    inds_cnfgrtns =  removeReferenceToConfigurations(intrsctn, ... %intersection structure version to edit
                                                     mask_pair, ... %index of the intersections in a pair of intersecting intersections
                                                     ind_cnddt_ln, ... %candiate lines which is being removed
                                                     inds_configurations_remove); %indice configuration to remove
                                                 
    try
        intrsctn.cnfgrtns{mask_pair}{ind_cnddt_ln} = inds_cnfgrtns;
    catch e
        rethrow(e);
        
    end
    
    if isempty(inds_cnfgrtns) 
       intrsctn = removeReferenceToCandidateLine(intrsctn, ... %intersection structure version to edit
                                         mask_pair, ... %index of the intersections in a pair of intersecting intersections
                                         ind_cnddt_ln... %candiate lines which is edited
                                         );
    end
end

function inds_cnfgrtns =  removeReferenceToConfigurations(intrsctn, ... %intersection structure version to edit
                                         mask_pair, ... %index of the intersections in a pair of intersecting intersections
                                         ind_cnddt_ln, ... %candiate lines which is edited
                                         inds_configurations_remove) %indice configuration to remove

    if isempty(ind_cnddt_ln)
        inds_cnfgrtns  = [];
        return;
    end
    
    try
        inds_cnfgrtns = intrsctn.cnfgrtns{mask_pair}{ind_cnddt_ln};
    catch e
        rethrow(e)
    end
    mask_cnfgrtn_kp = ~ismember(inds_cnfgrtns, inds_configurations_remove);
    inds_cnfgrtns = inds_cnfgrtns(mask_cnfgrtn_kp);
   
end

function intrsctn = removeReferenceToCandidateLine(intrsctn, ... %intersection structure version to edit
                                         mask_pair, ... %index of the intersections in a pair of intersecting intersections
                                         ind_cnddt_ln... %candiate lines which is edited
                                         )

   inds_cnddt_lns_keep = setdiff(1:length(intrsctn.cnddt_lns{mask_pair}), ind_cnddt_ln);
   if ~isempty(inds_cnddt_lns_keep)
        intrsctn.cnddt_lns{mask_pair} = intrsctn.cnddt_lns{mask_pair}(inds_cnddt_lns_keep);
        intrsctn.cnfgrtns{mask_pair} = intrsctn.cnfgrtns{mask_pair}(inds_cnddt_lns_keep);
   else
        intrsctn.cnddt_lns{mask_pair} = [];
   end
end

% function strokes_topology = reindexIntersections(ind_strk_pr,...                        
%                                       intersections,...
%                                       strokes_topology,...        
%                                       ind_intrsctn,...
%                                       intrsctn_updt,...
%                                       ind_cnddt)
%                           
%                           
%     ind_strk = intersections(ind_intrsctn).strokes_indices(ind_strk_pr);
%     
%     for j = 1:length(intrsctn_updt.cnddt_lns{ind_strk_pr})
%         candidate_line = strokes_topology(ind_strk).candidate_lines(intrsctn_updt.cnddt_lns{ind_strk_pr}(j));
%         
%         for jj = 1:length(intrsctn_updt.cnfgrtns{ind_strk_pr}{j})
%             ind_cnfgrtn = intrsctn_updt.cnfgrtns{ind_strk_pr}{j}(jj);
%             cnfgrtn = candidate_line.configurations(ind_cnfgrtn);
% 
%             mask_reindex = cnfgrtn.inds_intrsctns__mult_cnddts == ind_intrsctn;
% %             fprintf('Changes from %d to %d \n', cnfgrtn.inds_intrsctns__mult_cnddts_ind(mask_reindex), ind_cnddt);
%             cnfgrtn.inds_intrsctns__mult_cnddts_ind(mask_reindex) = ind_cnddt;
% 
%             candidate_line.configurations(ind_cnfgrtn) = cnfgrtn;
%         end
%         
%         strokes_topology(ind_strk).candidate_lines(intrsctn_updt.cnddt_lns{ind_strk_pr}(j)) = candidate_line;
%     end
% 
% end