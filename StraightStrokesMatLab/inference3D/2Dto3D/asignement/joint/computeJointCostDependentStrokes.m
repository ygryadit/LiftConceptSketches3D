function [max_cost,...
          strokes_topology,...
          num_joint_strokes,...
          list_joint_strokes] = ...
                computeJointCostDependentStrokes(cost_parent_stroke,...
                                      strks_visited,...
                                      strks_jnt_cnddt_lns,...
                                      strks_jnt_cnddt_lns_cnfgrtns,...
                                      configuration, ...
                                      strokes_topology, ...
                                      intersections)
     
    %% I. No dependent strokes:
     % We are here only if stroke 'cur_stroke_ind' is not a member of
     % strks_visited:
     if isempty(configuration.inds_jnts_strks)
        % There is no dependent strokes.
        max_cost = cost_parent_stroke;
        num_joint_strokes  = 1;
        return;
     end       

     %% II. All dependent strokes were already visited along the branch:
     
     % Find a subset of not accounted yet joint strokes:
     [inds_jnts_strks, inds_cnsdr] = setdiff(configuration.inds_jnts_strks, strks_visited);
     
     inds_intrsctns__mult_cnddts = configuration.inds_intrsctns__mult_cnddts(inds_cnsdr);
     inds_intrsctns__mult_cnddts_ind = configuration.inds_intrsctns__mult_cnddts_ind(inds_cnsdr);
     num_jnt_strks = length(inds_jnts_strks);     
     
     list_joint_strokes = inds_jnts_strks;

     
     if size(inds_jnts_strks,1) ~= 1
        inds_jnts_strks = inds_jnts_strks';
     end
      
     % All the dependent strokes are already taken into account:
     if num_jnt_strks == 0
        max_cost = cost_parent_stroke;
        num_joint_strokes  = 1;
        return;
     end       
     
     
     %% III. Accumulate costs from all dependent strokes:
     
     % For each intersection from the lines with multiple hypothesis:
%      fprintf('length(ind_strks_cnsdr) %d\n', length(ind_strks_cnsdr));
     max_cost_branch = cell(num_jnt_strks,1);
     num_elements_branch = cell(num_jnt_strks,1);
     
     for j = 1:num_jnt_strks 
        
        ind_intrsctng_strk  = inds_jnts_strks(j);
        ind_intrsctn        = inds_intrsctns__mult_cnddts(j);
        ind_intrsctn_vrsn   = inds_intrsctns__mult_cnddts_ind(j);
       
        if (ind_intrsctng_strk == 22) & (ind_intrsctn == 80) & (ind_intrsctn_vrsn == 1)
            
           disp(''); 
        end
%         
        intrsctn_cnddts3D = intersections(ind_intrsctn).cnddts3D(ind_intrsctn_vrsn);
        
        % Find which of the two intersecting strokes is not the current
        % stroke:
        mask_strk_cnfgrtn = intersections(ind_intrsctn).strokes_indices == ind_intrsctng_strk; 
            
        % Find which candidate geometry of intersecting stroke belongs
        % the intersection:        
        inds_cnddt_lns = intrsctn_cnddts3D.cnddt_lns{mask_strk_cnfgrtn};
        
        for g = 1:length(inds_cnddt_lns)
            ind_cnddt_g = inds_cnddt_lns(g);
            
            num_configs = length(intrsctn_cnddts3D.cnfgrtns{mask_strk_cnfgrtn}{g});
   
            for k = 1:num_configs
                
                ind_cnfgrtn = intrsctn_cnddts3D.cnfgrtns{mask_strk_cnfgrtn}{g}(k);

                configuration_brnch = strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_g).configurations(ind_cnfgrtn);
                
                if (isfield(configuration_brnch, 'pfj_assgnd') & (configuration_brnch.pfj_assgnd == strks_visited(1)))
                    %The branch was already visited and the joint cost does
                    %not need to be recomputed.                    
                    max_cost_ = configuration_brnch.p_full_joint;
                    num_joint_strokes_ = configuration_brnch.num_joint_strokes;
                    list_joint_strokes_ = configuration_brnch.list_jnt_strks;
                    
                elseif (configuration_brnch.p_full > 0.5)
%                     fprintf('Strk % d Line %d, configuration %d, cost full %.3f \n', ...
%                         ind_intrsctng_strk, ind_cnddt_g, ind_cnfgrtn, configuration_brnch.p_full);
                    
                    [max_cost_, ...
                     strokes_topology, ...
                     num_joint_strokes_, ...
                     list_joint_strokes_] = ...
                        computeJointCostDependentStrokes( [configuration_brnch.p_full],... 
                                                          [strks_visited inds_jnts_strks configuration.inds_jnts_strks],...
                                                          [strks_jnt_cnddt_lns          ind_cnddt_g],...
                                                          [strks_jnt_cnddt_lns_cnfgrtns ind_cnfgrtn],...
                                                          configuration_brnch, ...
                                                          strokes_topology, ...
                                                          intersections);
%                     fprintf('Done: Strk % d Line %d, configuration %d, cost full %.3f \n', ...
%                           ind_intrsctng_strk, ind_cnddt_g, ind_cnfgrtn, configuration_brnch.p_full);
                   
                else
                    max_cost_ = configuration_brnch.p_full; 
                    num_joint_strokes_ = 1; 
                    list_joint_strokes_ = configuration_brnch.list_jnt_strks;
                end
                
                max_cost_branch{j}(end+1) = max_cost_;
                num_elements_branch{j}(end+1) = num_joint_strokes_;
        
                strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_g).configurations(ind_cnfgrtn).p_full_joint = max_cost_;
                strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_g).configurations(ind_cnfgrtn).num_joint_strokes = num_joint_strokes_;             
                strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_g).configurations(ind_cnfgrtn).pfj_assgnd = strks_visited(1);  
                strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_g).configurations(ind_cnfgrtn).list_jnt_strks = list_joint_strokes_;
                
                list_joint_strokes = unique([list_joint_strokes, ...
                                          list_joint_strokes_]);
            end
    
            strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_g).max_cost = ...
                max(cat(1, strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_g).configurations(:).p_full_joint));
                

            strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_g).list_jnt_strks = unique(list_joint_strokes);    
        end
        
        
        costs_cndt_lns = cat(1, strokes_topology(ind_intrsctng_strk).candidate_lines(:).max_cost);
        
        strokes_topology(ind_intrsctng_strk).score = ...
                    max(costs_cndt_lns);
                
                
        [~, strokes_topology(ind_intrsctng_strk).score, ...
            strokes_topology(ind_intrsctng_strk).confidence] = ...
                assignCostAndConfidance(costs_cndt_lns);
                    
                    
        
            [max_cost_branch{j}, ind_max] = max(max_cost_branch{j});
            num_elements_branch{j} = num_elements_branch{j}(ind_max);        
     end
     
     num_joint_strokes = sum(cat(1,num_elements_branch{:})) + 1;
%      fprintf('cat(1,max_cost_branch{:}).*cat(1,num_elements_branch{:}) %.3f \n', cat(1,max_cost_branch{:}).*cat(1,num_elements_branch{:}))
     max_cost = sum([cost_parent_stroke; cat(1,max_cost_branch{:}).*cat(1,num_elements_branch{:})])./num_joint_strokes;

     
end
