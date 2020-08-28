function [candidate_lines] = findAllIntersectionsNearTheLine(cur_stroke, ...
                                            intersections,...
                                            candidate_lines,...
                                            cam_param,...
                                            depth_assigned,...
                                            strokes_topology)    
    %% Setup the constant threshold:
%     ratio_line_theshold = 1e-1;
    USE_3D_TOLERANCE = false;

    ratio_line_theshold = 0.1;

    %% Go over the candiate lines:
    num_candidate_lines = length(candidate_lines);
    
%     global last_added_stroke;
    
    for l_ind = 1:num_candidate_lines
        
%         try
%         mask_assigned_strokes_up_to_the_last_stroke = (cur_stroke.indcs_intrsctng_strks <= last_added_stroke) & ...
%                                                    depth_assigned(cur_stroke.indcs_intrsctng_strks)';
%         
%         ind_active_intersections = cur_stroke.indcs_intrsctns(mask_assigned_strokes_up_to_the_last_stroke)';  
%         catch e
%             rethrow(e);
%         end
        ind_active_intersections = cur_stroke.inds_intrsctns_eval_actv';
        
        intesetcion_coord3D = cat(1,intersections(ind_active_intersections).coordinates3D);
        
        if isempty(intesetcion_coord3D)
            continue;
        end
        
        point1_3Dline = candidate_lines(l_ind).origin - ...
                         0.5*candidate_lines(l_ind).dir;
                     
        point2_3Dline = candidate_lines(l_ind).origin + ...
                         0.5*candidate_lines(l_ind).dir;
        
        distances_eq = pointLineDistance3D(intesetcion_coord3D, ...
                                           point1_3Dline,...
                                           point2_3Dline);

        distances_proj = pointLineDistance3DPerspective(intesetcion_coord3D, ...
                                           point1_3Dline,...
                                           point2_3Dline,...
                                           cam_param);
            
        ratios_eq = distances_eq/candidate_lines(l_ind).length3D;
         
        if USE_3D_TOLERANCE
            tolerance = computeToleranceDistance(strokes_topology, candidate_lines(l_ind));
            ind_belong_line = find((distances_eq < tolerance) & (distances_proj < tolerance));
        else
                                       
            
            ratios_proj = distances_proj/candidate_lines(l_ind).length3D;

            ind_belong_line = find((ratios_eq < ratio_line_theshold) & (ratios_proj < ratio_line_theshold));
        end
        % Update the list of intersections belonging to the line (difference due to active intersections):
        [candidate_lines(l_ind).configurations(1).inds_intrsctns, ia, ic] = ...
            unique([ind_active_intersections(ind_belong_line), ...
                   candidate_lines(l_ind).configurations(1).inds_intrsctns], ...
                   'stable');
        
        % Update the list of active intersections belonging to the line:
        [candidate_lines(l_ind).configurations(1).inds_intrsctns__assigned, ~, ~] = ...
            unique([ind_active_intersections(ind_belong_line),...
                    candidate_lines(l_ind).configurations(1).inds_intrsctns__assigned],...
                    'stable');
        
        % Update the disatnce to intersections: not sure that this is used:
         candidate_lines(l_ind).configurations(1).p_intrsctns_dists = ...
             [1.0 - ratios_eq(ind_belong_line)', ...
              ones(size(candidate_lines(l_ind).configurations(1).p_intrsctns_dists))];
          
         candidate_lines(l_ind).configurations(1).p_intrsctns_dists = ...
             candidate_lines(l_ind).configurations(1).p_intrsctns_dists(ia);
         
%          fprintf('Intersections:'); disp(candidate_lines(l_ind).configurations(1).inds_intrsctns);
%          fprintf('Distances:'); disp(candidate_lines(l_ind).configurations(1).p_intrsctns_dists);
    end
end