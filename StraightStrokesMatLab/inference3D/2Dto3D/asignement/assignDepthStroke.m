function [strokes_topology, intersections] = ...
               assignDepthStroke(   strokes_topology,....
                                    intersections,...
                                    ind_cur_strk, ...                                    
                                    cam_param,...
                                    pairsInterInter,...
                                    do_assign_depth)
    
    %% Load global information:
    global fid;
    global DEBUG;
    
    global DISPLAY_INFO;
    global SHOW_FIGS;
    global filepath_sketch_img;
        
    global last_added_stroke;
    last_added_stroke = max(ind_cur_strk, last_added_stroke);
    
    fprintf(fid,'Current stroke = %d\n', ind_cur_strk);                                       

    %% Check input parameters:
    if ~exist('do_assign_depth', 'var')
        do_assign_depth = true;
    end 

   %% Indices of nodes of the current stroke:
    cur_stroke = strokes_topology(ind_cur_strk);
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
    
    
    if isempty(cur_stroke.inds_intrsctns_eval)
        return;
    end
     
    %% Figures:
    if SHOW_FIGS
        % Plot 3D: 


        fig_num_3D = 9;

        plotStrokesTopolgyIntersectionsTypes(...
                            strokes_topology,...
                            cur_stroke,...
                            intersections,...
                            fig_num_3D)
                        
       % Plot 2D: 
        tic                    
            plot2DCurStrokeIntersectingStokes(intersections, ...
                                              cur_stroke,...
                                              strokes_topology);
        ellapsed_time = toc;
        fprintf(fid,'\t\t Time to plot 2D intersections %.3f\n', ellapsed_time);
    end

    %% Try depth assignement:
    if ind_cur_strk > 128
        disp('');
    end
    try
        [strokes_topology, intersections] = ...
                    assignDepthStraightStroke(  cur_stroke,...
                                                intersections,...
                                                strokes_topology,...                                                                                     
                                                cam_param,...
                                                pairsInterInter,...
                                                false,...
                                                do_assign_depth);
    catch e
       if DEBUG
            fprintf(2,'Error in assigning depth to orthogonal line\n');
            rethrow(e)
       else
            rethrow(e)
       end
    end                           
                                 
    %% Figures:
   if SHOW_FIGS
        % Plot 3D: 
        fig_num_3D = 9;

        plotStrokesTopolgyIntersectionsTypes(...
                            strokes_topology,...
                            cur_stroke,...
                            intersections,...
                            fig_num_3D)
                        
       % Plot 2D: 
        tic    
        plot2DCurStrokeIntersectingStokes(intersections, ...
                                          cur_stroke,...
                                          strokes_topology);
        ellapsed_time = toc;
        fprintf(fid,'\t\t Time to plot 2D intersections %.3f\n', ellapsed_time);      
   end
   
   %% Additional infromation:
   
    if DISPLAY_INFO 
        if ~strokes_topology(ind_cur_strk).depth_assigned
            fprintf(fid, 'Stroke %d depth is not assigned\n', ind_cur_strk);
        else
            fprintf(fid, 'Stroke %d depth is assigned\n', ind_cur_strk);
        end
    end
    

end
