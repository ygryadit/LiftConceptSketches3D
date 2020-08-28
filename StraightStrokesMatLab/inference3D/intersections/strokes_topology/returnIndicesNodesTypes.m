% function [ind_intersections_with_prev_strokes,...
%           mask_indices_hypothesis,...
%           nodes_indices_active_others,...
%           mask_nodes_h_assigned_depth,...
%           indices_hypothesis_intersecting_strokes] = ...

% inds_intrsctns_eval
% inds_intrsctns_eval_actv
% inds_intrsctns_eval_mltpl_cnddts
%
% inds_intrsctng_strks_eval
% inds_intrsctng_strks_eval_actv
% inds_intrsctng_strks_eval_mltpl_cnddts


function [inds_intrsctns_eval,...  that are the intersections with previous strokes and assigned future strokes: inds_intrsctns_eval
          inds_intrsctns_eval_actv,...
          inds_intrsctns_eval_mltpl_cnddts,...
          inds_intrsctng_strks_eval,... that are the indices of the intersecting previous strokes and assigned future strokes: inds_intrsctng_strks_eval      
          inds_intrsctng_strks_eval_actv,...
          inds_intrsctng_strks_eval_mltpl_cnddts] = ...
                returnIndicesNodesTypes(cur_stroke, depth_assigned, intersections, UP_TO_LAST)

    global last_added_stroke;
    
    global USE_ONLY_LIKELY_INTERS;
    
    if USE_ONLY_LIKELY_INTERS
%         mask_likely = cat(1,intersections(cur_stroke.indcs_intrsctns).likely) | cat(1,intersections(cur_stroke.indcs_intrsctns).collinear);
        mask_likely = cat(1,intersections(cur_stroke.indcs_intrsctns).likely);
        
        cur_stroke.indcs_intrsctng_strks    = cur_stroke.indcs_intrsctng_strks( mask_likely );
        cur_stroke.indcs_intrsctns          = cur_stroke.indcs_intrsctns( mask_likely );
    end    
    
    if ~exist('UP_TO_LAST', 'var')
        UP_TO_LAST = false;
    end
    
    % Mask of intersections with previous strokes or up to last stroke:
    if ~UP_TO_LAST
        mask_intrsctns_prvs_strks = cur_stroke.indcs_intrsctng_strks <= cur_stroke.ind;
    else
        mask_intrsctns_prvs_strks = cur_stroke.indcs_intrsctng_strks <= last_added_stroke;
    end
    
    % Mask of intersections up to last assigned stroke:
    mask_assigned_strokes_up_to_the_last_stroke = (cur_stroke.indcs_intrsctng_strks <= last_added_stroke) & ...
                                                   depth_assigned(cur_stroke.indcs_intrsctng_strks);

    mask_set = mask_assigned_strokes_up_to_the_last_stroke | mask_intrsctns_prvs_strks;

    %Intersecting strokes:
    inds_intrsctng_strks_eval              = cur_stroke.indcs_intrsctng_strks(mask_set);  
    inds_intrsctng_strks_eval_actv         = inds_intrsctng_strks_eval(depth_assigned(inds_intrsctng_strks_eval));
    inds_intrsctng_strks_eval_mltpl_cnddts = inds_intrsctng_strks_eval(~depth_assigned(inds_intrsctng_strks_eval));
    
    %Intersections:
    inds_intrsctns_eval              = cur_stroke.indcs_intrsctns(mask_set);  
    inds_intrsctns_eval_actv         = inds_intrsctns_eval(depth_assigned(inds_intrsctng_strks_eval));
    inds_intrsctns_eval_mltpl_cnddts = inds_intrsctns_eval(~depth_assigned(inds_intrsctng_strks_eval));
    
    
  
%     indcs_strks_mltpl_cnddts = ...
%         cur_stroke.inds_intrsctng_strks_eval_mltpl_cnddts;
%     
%         %cur_stroke.indcs_intrsctng_strks(mask_intrsctns_prvs_strks);
% 
%     %% Commented out since any non resolved intersections with previous stroke should be taken into account:
% %     inds_intrsctns_eval    = inds_intrsctns_eval(depth_assigned(indcs_strks_mltpl_cnddts));
% %     mask_intrsctns_prvs_strks     = mask_intrsctns_prvs_strks & depth_assigned(cur_stroke.indcs_intrsctng_strks);
%     mask_indcs_intrsctns_prvs_strks_actv = depth_assigned(cur_stroke.inds_intrsctng_strks_eval);
% 
%     %%
%     % Other indices of the intersections that are active but not hypothesis:
%     indcs_intrsctns_actv_fllwng_strks    = cur_stroke.indcs_intrsctns( ...
%                                         depth_assigned ( ...
%                                             cur_stroke.indcs_intrsctng_strks ) );
%                                         
%     indcs_intrsctns_actv_fllwng_strks   = setdiff(indcs_intrsctns_actv_fllwng_strks, inds_intrsctns_eval);
    
%     disp('cur_stroke.indcs_intrsctns '); disp(cur_stroke.indcs_intrsctns);
%     disp('inds_intrsctns_eval: '); disp(inds_intrsctns_eval);
%     disp('mask_intrsctns_prvs_strks: '); disp(mask_intrsctns_prvs_strks);
%     disp('indcs_intrsctns_actv_fllwng_strks: '); disp(indcs_intrsctns_actv_fllwng_strks);
    
    
    
%     ind_nodes_hypothesis  = cur_stroke.inds_intrsctns_eval(depth_assigned(cur_stroke.indcs_intrsctng_strks(cur_stroke.mask_intrsctns_prvs_strks)));
%     ind_nodes_constraints = cur_stroke.indcs_intrsctns(depth_assigned(cur_stroke.indcs_intrsctng_strks));  
%     ind_strokes_constraints = cur_stroke.indcs_intrsctng_strks(depth_assigned(cur_stroke.indcs_intrsctng_strks));
%     
%     ind_nodes_active = find(intersections.mask_active);
%     ind_nodes_active = setdiff(ind_nodes_active, ind_nodes_constraints);    
%     
%      
%     
%     nodes_active.coordinates2D   = intersections.coordinates2D(ind_nodes_active, :);
%     nodes_active.coordinates3D   = intersections.coordinates3D(ind_nodes_active, :);
%     nodes_active.p_dist_str_segs = intersections.p_dist_str_segs(ind_nodes_active, :);
%     nodes_active.strokes_indices = intersections.strokes_indices(ind_nodes_active, :);
%     nodes_active.direction_types = intersections.direction_types(ind_nodes_active, :);
%     nodes_active.appear_number = mean(intersections.strokes_indices(ind_nodes_active, :),2);
%     nodes_active.ind_nodes_active = ind_nodes_active;
%     
%     
%     ind_non_coinciding_lines = nodes_active.direction_types;
%     % Only intersections of type 1:
%     nodes_active.coordinates2D = nodes_active.coordinates2D(ind_non_coinciding_lines,:);
%     nodes_active.coordinates3D = nodes_active.coordinates3D(ind_non_coinciding_lines,:);
%     nodes_active.p_dist_str_segs = nodes_active.p_dist_str_segs(ind_non_coinciding_lines,:);
%     nodes_active.strokes_indices = nodes_active.strokes_indices(ind_non_coinciding_lines,:);
%     nodes_active.direction_types = nodes_active.direction_types(ind_non_coinciding_lines, :);
%     nodes_active.appear_number = nodes_active.appear_number(ind_non_coinciding_lines, :);
%     nodes_active.ind_nodes_active = nodes_active.ind_nodes_active(ind_non_coinciding_lines);
%     
%     nodes_constraints.coordinates2D   = intersections.coordinates2D(ind_nodes_constraints, :);
%     nodes_constraints.coordinates3D   = intersections.coordinates3D(ind_nodes_constraints, :);
%     nodes_constraints.p_dist_str_segs = intersections.p_dist_str_segs(ind_nodes_constraints, :);
%     nodes_constraints.strokes_indices = intersections.strokes_indices(ind_nodes_constraints, :);
%     nodes_constraints.ind_nodes_constraints = ind_nodes_constraints; % among all the intersections
%     
%     nodes_constraints.appear_number = mean(intersections.strokes_indices(ind_nodes_constraints, :),2);
%     
%     nodes_constraints.mask_nodes_hypothesis = ismember(ind_nodes_constraints, ind_nodes_hypothesis); % among constraints
%     
%     nodes_constraints.strokes_indices = ind_strokes_constraints; % intersections.strokes_indices(ind_nodes_constraints,:);
    
%     nodes_constraints.length2D = intersections.length2D(ind_nodes_constraints,:);
%     nodes_constraints.threshold_upper_distance = intersections.threshold_upper_distance(ind_nodes_constraints,:);
end