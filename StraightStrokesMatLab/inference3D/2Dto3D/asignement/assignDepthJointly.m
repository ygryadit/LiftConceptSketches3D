function [strokes_topology, ...
          intersections] = assignDepthJointly(strokes_topology, ...
                                                 intersections,...
                                                 cam_param,...
                                                 candidate_line,...
                                                 cur_stroke,...
                                                 p_joint_configurations,...
                                                 max_cost,...
                                                 pairsInterInter)
    global SHOW_FIGS;
    %% Info:
    assigned_stroke = cur_stroke.ind;
%     if assigned_strokes == 37
%        disp('check'); 
%     end    

    %% I. Find the configuration that yeilds the maximum value:
    ind_mxmzng_cnfgrtn = find(p_joint_configurations == max_cost);
    
    % Assign the wining score to the curent stroke:
    try
        strokes_topology(cur_stroke.ind).score_alignment =  ...
            candidate_line.configurations(ind_mxmzng_cnfgrtn(1)).p_directional;
        strokes_topology(cur_stroke.ind).score_coverage = ...
            candidate_line.configurations(ind_mxmzng_cnfgrtn(1)).p_coverage;
        strokes_topology(cur_stroke.ind).planes_normals = ...
            candidate_line.configurations(ind_mxmzng_cnfgrtn(1)).planes_normals;
    catch e
        rethrow(e);
    end
      
    
    %% Add the intersection normal to each assigned stroke that intersects with the current stroke
    intrsctns__assigned = ...
        candidate_line.configurations(ind_mxmzng_cnfgrtn(1)).inds_intrsctns__assigned;
    for i = 1:length(intrsctns__assigned )
        ind_intrsctng_strk = setdiff(intersections(intrsctns__assigned(i)).strokes_indices, ...
                                     assigned_stroke);
        normal = cross(candidate_line.dir, strokes_topology(ind_intrsctng_strk).direction_vec);
        if norm(normal) > 0.5
            strokes_topology(ind_intrsctng_strk).planes_normals(end+1,:) = [normal ind_intrsctng_strk];
        end
    end
    
    %% II:
    % If there is one wining configuration, then find it, else consider all
    % the configurations:
    
    inds_configurations = getWinningConfigurations(p_joint_configurations, candidate_line);
    
    candidate_line_all_configurations = candidate_line;
    
    candidate_line.configurations = candidate_line.configurations(inds_configurations);
    
    
    %% III:
    % First find the set of intersections that are common for all the
    % configuartions of the optimal line:
    
    [intrsctns_cmmn, intrsctns_mlt_hpthss, strks_mlt_hpthss, inds_cnddt_intrsctn] = ...
        findCommonSubsetOfIntersections(candidate_line.configurations, intersections, strokes_topology);
    
    
    assign_table = table(intrsctns_mlt_hpthss, strks_mlt_hpthss, inds_cnddt_intrsctn);
    
    
    % All the potential intersections in the candiate line:
    ind_intersections_line = cat(2,candidate_line.configurations(:).inds_intrsctns);
        
    %ind_inter_nonactive = setdiff(cur_stroke.inds_intrsctns_prvs_strks, ind_intersections_line);
    %ind_inter_nonactive = setdiff(cur_stroke.inds_intrsctns_prvs_strks(cur_stroke.mask_indcs_intrsctns_prvs_strks_actv), ind_intersections_line);
    ind_inter_nonactive = setdiff(cur_stroke.inds_intrsctns_eval_actv, ind_intersections_line);
    
    ind_inter_active = intrsctns_cmmn; % these intersections we are sure should be assigned
    
    %There can be some non decided intersections.
 
    % Active intersections:
    vals_assign = num2cell(true(size(ind_inter_active)));
    [intersections(ind_inter_active).is_active] = vals_assign{:};
    
    % Non-active intersections:
    vals_assign = num2cell(false(size(ind_inter_nonactive)));
    [intersections(ind_inter_nonactive).is_active] = vals_assign{:};
    
    
    
    %% Clean up:
    [inds_strks_zero_cnddts, strokes_topology, intersections ] =  ...
                cleanUpTree(cur_stroke, ...
                            strokes_topology, ...
                            intersections,...
                            candidate_line_all_configurations.configurations);
    
    %% Assign depth:
    %Fill in the obtained 3D data for the newly assigned stroke:
    try
    [strokes_topology(cur_stroke.ind), intersections] = ...
        assignStrokeOneCandidateLine(strokes_topology(cur_stroke.ind),...
                                     intersections, ...
                                     cam_param, ...
                                     candidate_line);
        if isfield(strokes_topology(cur_stroke.ind), 'merged_with')
            %Define new variable:
            merged_with = setxor(strokes_topology(cur_stroke.ind).merged_with, cur_stroke.ind);
            
            
            for i = 1:length(merged_with)
               ind_s =  merged_with(i);
               strokes_topology(ind_s).indcs_intrsctns = [];
                [strokes_topology(ind_s), intersections] = ...
                    assignStrokeOneCandidateLine(strokes_topology(ind_s),...
                                         intersections, ...
                                         cam_param, ...
                                         candidate_line);

            end
            disp('Assigned depth addtiional strokes');
        end
    catch e
       rethrow(e);
    end
    global last_added_stroke;
    strokes_topology(cur_stroke.ind).assigned = last_added_stroke;                                    
    
    
      if SHOW_FIGS
        fig_num_3D = 9;
        
        plotStrokesTopolgyIntersectionsTypes(...
                            strokes_topology,...
                            cur_stroke,...
                            intersections,...
                            fig_num_3D) 
        
        plot2DCurStrokeIntersectingStokes(intersections, ...
                                          cur_stroke,...
                                          strokes_topology);
      
      end
      
    %% Clean
    for j = 1:length(intrsctns_mlt_hpthss)
        
            if isempty(strokes_topology(assign_table.strks_mlt_hpthss(j)).candidate_lines)
                continue;
            end
            [strokes_topology,...
               intersections] = assignDepthJointStroke( intersections,...
                                                        strokes_topology,...
                                                        assign_table,...
                                                        j,...
                                                        cam_param,...
                                                        pairsInterInter);   
             if ~isempty(assign_table.intrsctns_mlt_hpthss(j))
                 intersections(assign_table.intrsctns_mlt_hpthss(j)).cnddts3D = []; %clean 
             end    
        
    end
    
    %% Create candidate lines and try assign strokes that previously had no candidate lines
    global folder_save;
    
    for ind_strk_zero_cnddts = inds_strks_zero_cnddts
          if (~isfield(strokes_topology(ind_strk_zero_cnddts), 'candidate_lines') || ...
                isempty(strokes_topology(ind_strk_zero_cnddts).candidate_lines)) && ...
                (~strokes_topology(ind_strk_zero_cnddts).depth_assigned) && ...
                (ind_strk_zero_cnddts ~= last_added_stroke)
            
              [strokes_topology, intersections] = ...
                   assignDepthStroke(strokes_topology,....
                                     intersections,...
                                     ind_strk_zero_cnddts, ...
                                     cam_param,...
                                     pairsInterInter,...
                                     true);
          end
    end
    
    
    
%     for j = 1:length(cur_stroke.inds_intrsctns_eval_mltpl_cnddts)
%        ind_intrsctn = cur_stroke.inds_intrsctns_eval_mltpl_cnddts(j);
%        mask_ismeber = ismember(intrsctns_mlt_hpthss, ind_intrsctn);
%        if mask_ismeber
%           % Assign
%           [strokes_topology,...
%            intersections] = assignDepthJointStroke( intersections,...
%                                                     strokes_topology,...
%                                                     assign_table,...
%                                                     find(mask_ismeber),...
%                                                     cam_param);
%        else
%            % Update stuctures
%            [strokes_topology, ...
%            intersections] = removeNonFeasibleIntersectionPositions(...
%                 strokes_topology, ...
%                 intersections,...
%                 candidate_line_all_configurations,... % the optimal candidate line (with all configurations)
%                 cur_stroke,...
%                 j); 
%        end
%     end
    

%     for i = 1:length(intrsctns_mlt_hpthss)
%         
%         ind_intrsctng_stk = find(intersections(intrsctns_mlt_hpthss(i)).strokes_indices == strks_mlt_hpthss(i));
%     
%         if (strokes_topology(strks_mlt_hpthss(i)).depth_assigned)
%             continue;
%         end
%         
%         inds_cnddt_ln = intersections(intrsctns_mlt_hpthss(i)).cnddts3D(inds_cnddt_intrsctn(i)).cnddt_lns{ind_intrsctng_stk};
%         
%         ind_configurations_with_activated_intersection = cell2mat(intersections(intrsctns_mlt_hpthss(i)).cnddts3D(inds_cnddt_intrsctn(i)).cnfgrtns{ind_intrsctng_stk});
%         % [Todo] desactivate some intersections if there are such in excluded
%         % configurations!
%         
%         % Keep only configuratins with this stroke:
%         %%% Generally, they can be more than one inds_cnddt_ln
%         if length(inds_cnddt_ln) > 1            
%             warning('maybe here we should do something else I am not sure');
%             continue;
%         end
%         
% %         strokes_topology(strks_mlt_hpthss(i)).candidate_lines(inds_cnddt_ln).configurations = ...
% %             strokes_topology(strks_mlt_hpthss(i)).candidate_lines(inds_cnddt_ln).configurations(ind_configurations_with_activated_intersection);
%   
%         candidate_line_joint = strokes_topology(strks_mlt_hpthss(i)).candidate_lines(inds_cnddt_ln);
%         candidate_line_joint.configurations = candidate_line_joint.configurations(ind_configurations_with_activated_intersection);
%         
%         stroke_joint = strokes_topology(strks_mlt_hpthss(i));
%         stroke_joint.ind = strks_mlt_hpthss(i);
%         
%         [stroke_joint.inds_intrsctns_eval,...
%          stroke_joint.inds_intrsctns_eval_actv,...
%          stroke_joint.inds_intrsctns_eval_mltpl_cnddts,...
%          stroke_joint.inds_intrsctng_strks_eval,...
%          stroke_joint.inds_intrsctng_strks_eval_actv,...
%          stroke_joint.inds_intrsctng_strks_eval_mltpl_cnddts] = ...        
%             returnIndicesNodesTypes(stroke_joint, ...
%                             cat(1, strokes_topology(:).depth_assigned));
%         
%         
%         %temporary:
%         p_joint_configurations_ = cat(1,candidate_line_joint.configurations(:).p_full_joint);
%         max_cost_joint = max(cat(1,candidate_line_joint.configurations(:).p_full_joint));
%          
%         [strokes_topology,...
%          intersections,...
%          vals] =...
%                 assignDepthJointly(strokes_topology, ...
%                                    intersections, ...
%                                    cam_param, ...
%                                    img, ...
%                                    candidate_line_joint, ...
%                                    strokes_topology(strks_mlt_hpthss(i)).candidate_lines,...
%                                    stroke_joint, ...
%                                    p_joint_configurations_,...
%                                    max_cost_joint); 
%         assigned_strokes((end+1):(end+length(vals))) = vals;
%         
%         if ~isempty(intrsctns_mlt_hpthss(i))
%             intersections(intrsctns_mlt_hpthss(i)).cnddts3D = []; %clean 
%         end
%     end
%     
    
    

        
  
end