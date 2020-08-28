function [strokes_topology, intersections] = ...
                assignDepthStraightStroke(  cur_stroke,...
                                    intersections,...
                                    strokes_topology,...                                                                                     
                                    cam_param,...
                                    pairsInterInter,...
                                    ALREADY_VISITED,...
                                    do_assign_depth)
    %% Globals
    global DISPLAY_INFO;      
    global SHOW_FIGS;
    global DO_PRUNE;
    
    global max_cost_threshold;
    global confidence_threshold;
    
    %% Check input variables:
    if ~exist('ALREADY_VISITED', 'var')
        ALREADY_VISITED = false;
    end
    
    if ~exist('do_assign_depth', 'var')
        do_assign_depth = true;
    end
        
    %% Directional prior:
    
    if ismember(cur_stroke.line_group, 1:3)
        % Direction towards vanishing lines:    
        direction_prior = getDirectionVec(cur_stroke.line_group);
    else
        direction_prior = [];
    end
   

    
    %% Find all the candidate lines:
    if ~ALREADY_VISITED
      expectedNumberCandidateLines = computeExpectedNumberCandidateLines(cur_stroke, strokes_topology, intersections);
      
      if DO_PRUNE
          [strokes_topology,...
           intersections] = pruneSolutionIfRequiredSimplified(cur_stroke,...
                                                      strokes_topology,...
                                                      intersections,...
                                                      pairsInterInter,...
                                                      cam_param,...
                                                      expectedNumberCandidateLines);
 
      end
     ind_cur_strk = cur_stroke.ind;
     cur_stroke = strokes_topology(cur_stroke.ind);
     cur_stroke.ind = ind_cur_strk;
   
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
                                    
                                    
     if SHOW_FIGS
        % Plot 3D: 
        fig_num_3D = 9;

        plotStrokesTopolgyIntersectionsTypes(...
                            strokes_topology,...
                            cur_stroke,...
                            intersections,...
                            fig_num_3D)
                        
       % Plot 2D: 

        plot2DCurStrokeIntersectingStokes(intersections, ...
                                          cur_stroke,...
                                          strokes_topology);


     end
   
      candidate_lines = ...
                findCandidateLines(  strokes_topology,...
                                     intersections,...
                                     pairsInterInter,...
                                     cur_stroke,...
                                     direction_prior,...
                                     cam_param);
    
       if isempty(candidate_lines)
            strokes_topology(cur_stroke.ind).depth_assigned = false;
            vals_assign = num2cell(false(size(cur_stroke.indcs_intrsctns)));
            [intersections(cur_stroke.indcs_intrsctns).is_active] = vals_assign{:};        
            return;
        end
      if cur_stroke.ind == 71
            disp('check');
      end
       
%       if DO_PRUNE
%           [candidate_lines,...
%               strokes_topology,...
%               intersections] = pruneSolutionIfRequired(cur_stroke,...
%                                                       candidate_lines,...
%                                                       strokes_topology,...
%                                                       intersections,...
%                                                       pairsInterInter,...
%                                                       cam_param);
%  
%       end
      
%         while ((length(candidate_lines) > thr_max_num_lines) || ((num_configurations) > thr_max_num_config)  ) && ...
%                 (confidence_threshold_reached > confidence_threshold_min) &&  ~costs_low
%             
%             confidence_threshold = confidence_step*confidence_threshold; 
%             
%             [candidate_lines,...
%                intersections,...
%             strokes_topology,...
%             confidence_threshold_reached,...
%             costs_low]...
%                 = assignTheOldestStrokeFirst(cur_stroke.ind,...
%                                         intersections,...
%                                         strokes_topology,...                                                                                     
%                                         cam_param,...
%                                         img,...
%                                         pairsInterInter);
%              confidence_threshold = confidence_threshold_old;                    
%             disp('Thresholds reduction');
%         end
        

        %% Add intersections into configurations of intersecting strokes:
        
        

        % Creates new configurations for new strokes, includign the
        % intersrction with the current stroke. It fills in the
        % intersections structure, including only references to the
        % intersecting strokes.
        [strokes_topology, intersections] = addIntersectionsIntoJointStrokes(  cur_stroke,...
                                                        candidate_lines,...
                                                        strokes_topology,...
                                                        intersections,...
                                                        cam_param,...
                                                        pairsInterInter);
                                                    
                                                    
        % Update directional priors and normals for the intersecting
        % strokes in the newly created configurations:
        try
            strokes_topology = updateDirectionalPriorsNormals( ...
                                                        strokes_topology,...
                                                        intersections,...
                                                        candidate_lines,...
                                                        cur_stroke.ind);
        catch e
           rethrow(e);
        end
        global last_added_stroke;
        strokes_topology(cur_stroke.ind).created = last_added_stroke; 
        
    else
        candidate_lines = strokes_topology(cur_stroke.ind).candidate_lines;
    end
    
    
     if isempty(candidate_lines)
         return;
     end
    
     for i = 1:length(candidate_lines)
        candidate_lines(i)  = computeNormalIntersectionsConfigurations(candidate_lines(i),...
                                                strokes_topology,...
                                                intersections,...
                                                cur_stroke);
     end                                       
    %% Line reliability based on the angle:
    try
        pDirectional = computeDirectionalCosts(candidate_lines,...
                                                direction_prior,...
                                                cur_stroke,...
                                                pairsInterInter,...
                                                strokes_topology,...
                                                intersections);
    catch e
        rethrow(e);
    end
    
    %% Compute the rest of cost as well as joint costs:
    tic
    [candidate_lines, p_full_cnddt_lns, p_joint_cnddt_lns, p_joint, strokes_topology] = ...
                                              computeJointCostsV2(  strokes_topology,...
                                                                    intersections,...
                                                                    candidate_lines,...
                                                                    pDirectional,...
                                                                    cur_stroke);
    elapsed_time = toc;
    fprintf('Time to compute joint cost = %.3f\n', elapsed_time);
    %% Sort the candidate lines:
    if ~ALREADY_VISITED
        [p_joint_cnddt_lns, inds_sorted] = sort(p_joint_cnddt_lns, 'descend');
        candidate_lines = candidate_lines(inds_sorted);
        p_joint = p_joint(inds_sorted);
    end
    %% Plot:

    if SHOW_FIGS
        fig_num = 11;
        num_plot = min(length(candidate_lines), 5);
        plotCandidateLines(candidate_lines(1:num_plot), ...
                           strokes_topology, ...
                           cur_stroke, ...
                           intersections, ...
                           fig_num);
    end

    
    %% Select the optimal line: 
    
    [ind_line_most_probable, max_cost, strokes_topology(cur_stroke.ind).confidence] = ...
        assignCostAndConfidance(p_joint_cnddt_lns);
    
    
    strokes_topology(cur_stroke.ind).score = max_cost;
    %Cost of the current stroke solely:
    max_cost_stroke = p_full_cnddt_lns(ind_line_most_probable);
     
    global DELAY_ASSIGN
    
    strokes_topology(cur_stroke.ind).num_candidate_lines_all = length(candidate_lines);
    
    if DELAY_ASSIGN 
%         if ~ALREADY_VISITED
            if doNotAssignDepthValueToLine(max_cost,...
                  strokes_topology(cur_stroke.ind).confidence,...
                  max_cost_stroke,...
                  length(p_joint_cnddt_lns),...
                  candidate_lines) || ~do_assign_depth

                  if ALREADY_VISITED
                      fprintf('Not assigned %d: max_cost = %.3f, confidence = %.3f\n',...
                          cur_stroke.ind, max_cost, strokes_topology(cur_stroke.ind).confidence)
                      
                      strokes_topology(cur_stroke.ind).candidate_lines = candidate_lines;
                      strokes_topology(cur_stroke.ind).score    = max_cost;
                      return;
                  end


                    max_cost = min(max_cost, max_cost_stroke);
                    strokes_topology(cur_stroke.ind).candidate_lines = candidate_lines;
                    strokes_topology(cur_stroke.ind).score    = max_cost;

                    

                    intersections = assignCandidateIntersectionSecondStroke(candidate_lines, intersections, cur_stroke.ind);
                    
                    strokes_topology(cur_stroke.ind).num_candidate_lines = length(strokes_topology(cur_stroke.ind).candidate_lines);
                                            
                    global USE_WEAK_REJECT;
                    if USE_WEAK_REJECT
                        strokes_topology(cur_stroke.ind).num_candidate_lines_before_trim = length(strokes_topology(cur_stroke.ind).candidate_lines);


                        if max_cost > max_cost_threshold
                            candidate_lines_keep = find(p_joint_cnddt_lns > max_cost*0.8);
                            candidate_lines_remove = find(p_joint_cnddt_lns <= max_cost*0.8);
                        else
                            if length(p_joint_cnddt_lns) < 10
                               candidate_lines_keep = 1:length(p_joint_cnddt_lns);
                               candidate_lines_remove = [];
                            else    
                               candidate_lines_keep = find(p_joint_cnddt_lns > 0.2);
                               candidate_lines_remove = find(p_joint_cnddt_lns <= 0.2);
                            end
                        end

                        if SHOW_FIGS
                            num_plot = min(length(candidate_lines_keep), 5);
                            plotCandidateLines(candidate_lines(candidate_lines_keep(1:num_plot)), ...
                                               strokes_topology, ...
                                               cur_stroke,...
                                               intersections, 51);
                        end

                        if ~isempty(candidate_lines_remove)
                            [strokes_topology, intersections] = ...
                                removeCandidateLine(strokes_topology, ...
                                                    cur_stroke,...
                                                    candidate_lines_remove,...
                                                    intersections);
                        end

                        strokes_topology(cur_stroke.ind).num_candidate_lines = length(strokes_topology(cur_stroke.ind).candidate_lines);
                        strokes_topology(cur_stroke.ind).num_candidate_lines_after_trim = strokes_topology(cur_stroke.ind).num_candidate_lines;
                    end
                    
                    if DISPLAY_INFO
                        fprintf('Depth is not assigned,\n Numer of candidate lines = %d,\n', ...
                                strokes_topology(cur_stroke.ind).num_candidate_lines);
                        fprintf('Best score %.3f, confidence = %.3f\n',...
                                max_cost,...
                                strokes_topology(cur_stroke.ind).confidence);
                    end
                    return;
            end
            
            %% Select the stroke with closest in time intersections.
            if strokes_topology(cur_stroke.ind).confidence < confidence_threshold
               % Find the strokes those intersections are with strokes closer in
               % time.

               inds = find(p_joint_cnddt_lns > 0.98);
               if ~isempty(inds )
                   vals_dist = zeros(length(inds),1);
                   for i = 1:length(inds)
                      ind_config = find([candidate_lines(inds(i)).configurations(:).p_full_joint] == ...
                                        candidate_lines(inds(i)).max_cost);
                      inds_intrsctns__ = candidate_lines(inds(i)).configurations(ind_config).inds_intrsctns;
                      strokes_indices = [intersections(inds_intrsctns__).strokes_indices];
                      strokes_indices = setdiff(strokes_indices(:), cur_stroke.ind);
                      vals_dist(i) = mean(abs(strokes_indices)-cur_stroke.ind);
                   end
                   inds_line_most_probable = inds(vals_dist==min(vals_dist));
                   max_costs = candidate_lines(inds_line_most_probable).max_cost;
                   [max_cost, ind] = max(max_costs);
                   ind_line_most_probable = inds_line_most_probable(ind);


                   ind_line_most_probable = ind_line_most_probable(1);
                   max_cost = max_cost(1);
               end
            end
%         end
    end
    
    %% Assign depth value and remove axilary edges from the graph.
    
    

    
    
    fprintf('Assigned %d: max_cost = %.3f, confidence = %.3f\n',...
                          cur_stroke.ind, max_cost, strokes_topology(cur_stroke.ind).confidence)
    candidate_line = candidate_lines(ind_line_most_probable);
        
    [strokes_topology, intersections] = assignDepthJointly(strokes_topology,...
                                                               intersections,...
                                                               cam_param,...
                                                               candidate_line,...
                                                               cur_stroke,...
                                                               p_joint{ind_line_most_probable},...
                                                               max_cost,...
                                                               pairsInterInter);    
end






                                    
                                    