function  [strokes_topology, ...
           intersections,...
           strks_cleaned] = removeNonFeasibleIntersectionPositions(...
                strokes_topology, ...
                intersections,...
                all_configurations,... % the optimal candidate line 
                cur_stroke,...
                ind_intrsctn_update,...
                ind_intrsctng_strk) 
            

strks_cleaned = ind_intrsctng_strk;           
% assign_struct.intrsctns_mlt_hpthss
% assign_struct.strks_mlt_hpthss
% assign_struct.inds_cnddt_intrsctn            
            

% Iterate over all the intersections with multiple hypothesis:
% for j = 1:length(cur_stroke.inds_intrsctns_eval_mltpl_cnddts)
    %Find all the versions that can be it the current selected candidate
    %line:
    
    intrsctn_vrsn_kp = [];
    inds_confgrtns_with_cnddts_keep = [];
%     ind_intrsctn_update %= cur_stroke.inds_intrsctns_eval_mltpl_cnddts(j)
%     ind_intrsctng_strk% =  cur_stroke.inds_intrsctng_strks_eval_mltpl_cnddts(j)
    
    for ci = 1:length(all_configurations)
        mask_ismember = ismember(all_configurations(ci).inds_intrsctns__mult_cnddts, ind_intrsctn_update);
        if sum(mask_ismember)
            intrsctn_vrsn_kp(end+1) = all_configurations(ci).inds_intrsctns__mult_cnddts_ind(mask_ismember);
            inds_confgrtns_with_cnddts_keep(end+1) = ci;
        end
    end
    
    % for i = 1:length(inds_confgrtns_with_cnddts_keep)
        % if isempty(all_configurations(inds_confgrtns_with_cnddts_keep(i)).p_full_joint)
            % all_configurations(inds_confgrtns_with_cnddts_keep(i)).p_full_joint = ...
                % all_configurations(inds_confgrtns_with_cnddts_keep(i)).p_full;
        % end
    % end
    if length(intrsctn_vrsn_kp) > 1
        costs = cat(1, all_configurations(inds_confgrtns_with_cnddts_keep).p_full_joint);
        [~,ind_max] = max(costs);
        intrsctn_vrsn_kp_max = intrsctn_vrsn_kp(ind_max);
        intrsctn_vrsn_kp_max = intrsctn_vrsn_kp_max(1);
    elseif length(intrsctn_vrsn_kp) == 1
        intrsctn_vrsn_kp_max = intrsctn_vrsn_kp; 
    end
    
    intrsctn_vrsn_kp = unique(intrsctn_vrsn_kp);
    
    if ~isempty(intrsctn_vrsn_kp)
        intersections(ind_intrsctn_update).coordinates3D = intersections(ind_intrsctn_update).cnddts3D(intrsctn_vrsn_kp_max).coordinates3D;
    else 
        intersections(ind_intrsctn_update).is_active = false; 
    end
    
    
     %% Change the remaining option to be the intersection with assigned stroke but possibly no decision
    all_versions = 1:length(intersections(ind_intrsctn_update).cnddts3D);
    inds_cnddts_remove = setdiff(all_versions, intrsctn_vrsn_kp);
    
    ind_intrsctng_strk_pair = find(intersections(ind_intrsctn_update).strokes_indices == ind_intrsctng_strk);
    
%     %1. No versions of the intersection to keep.
%     if isempty(intrsctn_vrsn_kp)
%        %Is that even feasible? 
%        intersections(ind_intrsctn_update).cnddts3D = [];
%        intersections(ind_intrsctn_update).is_active = false;
%        return;
%     end
        
    %2. There is only one version of the intersection to keep:
    inds_cnd_lns_update_all = [];
    for ii = 1:length(intrsctn_vrsn_kp)
        
        inds_cnd_lns_update = intersections(ind_intrsctn_update).cnddts3D(intrsctn_vrsn_kp(ii)).cnddt_lns{ind_intrsctng_strk_pair};
        inds_cnd_lns_update_all = [inds_cnd_lns_update_all , inds_cnd_lns_update];
        for iclj = 1:length(inds_cnd_lns_update)
         
            inds_configurations_update = intersections(ind_intrsctn_update).cnddts3D(intrsctn_vrsn_kp(ii)).cnfgrtns{ind_intrsctng_strk_pair}{iclj};
            if isempty(inds_configurations_update)
                continue;
            end
            try
                
            [strokes_topology,...
             intersections] = updateConfigurations( ...
                                  ind_intrsctng_strk,...
                                  inds_cnd_lns_update(iclj),...
                                  inds_configurations_update,...
                                  ind_intrsctn_update,...
                                  strokes_topology,...
                                  intersections);
            catch e
                rethrow(e);
            end
          intersections(ind_intrsctn_update).cnddts3D(intrsctn_vrsn_kp(ii)).cnfgrtns{ind_intrsctng_strk_pair}{iclj} = 1:length(inds_configurations_update);
%           intersections(ind_intrsctn_update).coordinates3D = intersections(ind_intrsctn_update).cnddts3D(intrsctn_vrsn_kp).coordinates3D;  
%           intersections(ind_intrsctn_update).cnddts3D = [];
        end
        
    end
    
    
    
    
    
    %% Remove all other vesions and associted configurations of intersecting strokes
    % here the configurations are just removed
%     all_versions = 1:length(intersections(ind_intrsctn_update).cnddts3D);
%     inds_cnddts_remove = setdiff(all_versions, intrsctn_vrsn_kp);
%     
%     ind_intrsctng_strk_pair = find(intersections(ind_intrsctn_update).strokes_indices == ind_intrsctng_strk);
    
    if intersections(ind_intrsctn_update).is_active == 1
        % Remove all the candidate lines that did not have configuration with the assigned intersection:
        all_the_candidate_lines = 1:length(strokes_topology(ind_intrsctng_strk).candidate_lines);
        inds_cnddt_lns_clean = setdiff(all_the_candidate_lines, inds_cnd_lns_update_all);
        
        for k = inds_cnddt_lns_clean
%             try
              if isempty(strokes_topology(ind_intrsctng_strk).candidate_lines)
                  continue;
              end
              inds_configurations_remove = 1:length(strokes_topology(ind_intrsctng_strk).candidate_lines(k).configurations);
%             catch e
%                 rethrow(e);
%             end
           [strokes_topology,...
                intersections, ...
                strks_cleaned_] = removeConfigurations(ind_intrsctng_strk,...
                                                      k,...
                                                      inds_configurations_remove,...
                                                      strokes_topology,...
                                                      intersections,...
                                                      cur_stroke.ind);
                strks_cleaned = unique([strks_cleaned,strks_cleaned_]);
                
                
        end
%          [strokes_topology,...
%                 intersections, ...
%                 strks_cleaned] = cleanCandidateLines(strokes_topology,...
%                              intersections,...
%                              cur_stroke.ind,...   
%                              ind_intrsctng_strk,...
%                              inds_cnddt_lns_clean);
                         
    else
        for icr = inds_cnddts_remove
            try
                inds_cnd_lns = intersections(ind_intrsctn_update).cnddts3D(icr).cnddt_lns{ind_intrsctng_strk_pair};
            catch e
                rethrow(e)
            end
            %Iterate over candidate lines:
            for iclj = 1:length(inds_cnd_lns)
%                 if isempty(intersections(ind_intrsctn_update).cnddts3D(icr).cnfgrtns{ind_intrsctng_strk_pair})
%     %                warning('The empty candidate lines should be cleaned instead of ignoring them.')
%     %                continue; 
%                      error('');
%                 end

                inds_configurations_remove = intersections(ind_intrsctn_update).cnddts3D(icr).cnfgrtns{ind_intrsctng_strk_pair}{iclj};
                
                if isempty(inds_configurations_remove)
                   continue; 
                end
                
                try
                    [strokes_topology,...
                        intersections,...
                        strks_cleaned_] = removeConfigurations(ind_intrsctng_strk,...
                                                              inds_cnd_lns(iclj),...
                                                              inds_configurations_remove,...
                                                              strokes_topology,...
                                                              intersections,...
                                                              cur_stroke.ind);
                         strks_cleaned = unique([strks_cleaned,strks_cleaned_]);
                catch e
                    rethrow(e);
                end
    
            end
        end
    
    %By this point all the non-possible configurations from the
    %intersecting strokes are removed and indixing is correct.
    end
    
   
    
   
   
    
    intersections(ind_intrsctn_update).cnddts3D = []; %clean 
   
    %3. There are several versions of the intersection to keep:

    % Remove the intersection version:             
%     intersections(ind_intrsctn_update).cnddts3D = ...
%         intersections(ind_intrsctn_update).cnddts3D(intrsctn_vrsn_kp);


end








