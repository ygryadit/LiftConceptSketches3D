function [candidate_lines,...
          strokes_topology,...
          intersections] = ...
                pruneSolutionIfRequired(cur_stroke,...
                                        candidate_lines, ...
                                        strokes_topology, ...
                                        intersections,...
                                        pairsInterInter,...
                                        cam_param)
    % Load global parameters:
    global confidence_threshold;
    global thr_max_num_lines;
    global confidence_threshold_min;
    global max_cost_threshold;
        
    %Intialise:
    confidence_threshold_reached = confidence_threshold;
    costs_low = false;
    
    %Save:
    confidence_threshold_old = confidence_threshold;

    %Compute number of configurations for the current stroke:
    num_configurations = 0;
    for iii = 1:length(candidate_lines)
        num_configurations = num_configurations + length(candidate_lines(iii).configurations);
    end
    

    %Find the set of all dependent strokes:
    try
    [inds_non_assgnd_dpndnt_strks] = findNonAssgndDpndntStrks(cur_stroke,...
                                                 candidate_lines,...
                                                 strokes_topology,...
                                                 intersections)
    catch
       disp(''); 
    end
    % Prun till the criteria are sutisfied:    
    while ((length(candidate_lines) > thr_max_num_lines) || ...
            (checkNumConfigurationsAnyStroke(strokes_topology, candidate_lines, inds_non_assgnd_dpndnt_strks))  ) && ...
            (confidence_threshold_reached > confidence_threshold_min) 
        
        disp('Thresholds reduction');
        
        try
            [inds_non_assigned_strks,...
            max_costs,...
            confidence_vals,...
            confidence_threshold_reached] = ...
                    getNonAssigenedStrokesData(strokes_topology, max_cost_threshold, inds_non_assgnd_dpndnt_strks);
        catch e
            confidence_threshold = confidence_threshold_old;   
            rethrow(e);
        end
        
        fprintf('inds_non_assigned_strks:'); disp(inds_non_assigned_strks);
        
         if isempty(inds_non_assigned_strks) | ...
           (confidence_threshold_reached < confidence_threshold_min)  | ...
           (strokes_topology(inds_non_assigned_strks(1)).score < max_cost_threshold)
            break;
        end

       [strokes_topology,...
        intersections,...
        candidate_lines] = ...
                assignOnePrevStrokeCheckCandidateLines(cur_stroke.ind,...
                    strokes_topology,...
                    intersections,...
                    pairsInterInter,...
                    cam_param,...
                    inds_non_assigned_strks,...
                    max_costs,...
                    confidence_vals);
    end

    confidence_threshold = confidence_threshold_old;   

    end