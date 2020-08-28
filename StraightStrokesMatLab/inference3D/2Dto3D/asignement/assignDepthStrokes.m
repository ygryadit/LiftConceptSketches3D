% Assign depth to one stroke in a time.
% Revisites previously assigned strokes.


function [strokes_topology,....
          intersections] = ...
            assignDepthStrokes(strokes_topology,....
                              intersections,...
                              ind_first_stroke,...
                              folder_save,...
                              cam_param,...
                              pairsInterInter)
    %% =============== Initialise filepaths ===============================       
    
    global SAVE_GIFS
    global fid;
    global DEBUG;
    SAVE_GIFS = false;
    
    num_strokes = length(strokes_topology);    
    
    %% =============== Compute depth ======================================
    for ind_assign = (ind_first_stroke+1):num_strokes

        %% Skip non-straight strokes:
        if strokes_topology(ind_assign).primitive_type ~= 0
            % Skip for now non-line primitives.
            continue;
        end
       
        if ind_assign == 46
            disp('check');
        end
        
        %% Fill in the depth information for one more stroke:
        tic
        try
            [strokes_topology, intersections] = ...
                          assignDepthStroke(strokes_topology,....
                                           intersections,...
                                           ind_assign,...
                                           cam_param,...
                                           pairsInterInter);
        catch e
           if DEBUG
                fprintf(2,'Error in assigning depth to orthogonal line\n')  
                rethrow(e)
           else
                rethrow(e)
           end
        end
        elapsedTime = toc;
        fprintf(fid, 'Stroke %d / %d reconstruction: %.3f seconds\n', ...
                     ind_assign, num_strokes, elapsedTime);
        
        %% Save the current reconstruction result:
        try
            saveDrawingAsOBJ(strokes_topology, intersections, folder_save, 'confident');
            saveDrawingAsOBJSingleObject(strokes_topology, folder_save, 'confident');
            saveJSONReconstruction(strokes_topology, intersections, cam_param, folder_save, 'confident');
        catch e
           if DEBUG
                fprintf(2,'Error in saving data\n')     
                rethrow(e)
           else
                rethrow(e)
           end
        end
        
        %% Update the likelyhood of the intersections:
%           try                                         
%                 intersections = updateListLikelyIntersections(...
%                                 pairsInterInter,...
%                                 intersections,...
%                                 pairsInterAttrac,...
%                                 points_attraction_to_pair,...
%                                 collinear_strokes_pairs,...
%                                 img);   
%           catch
%              fprintf(2, 'Error updating likely nodes\n');
%           end
        
        %% ============= Save GIFs ========================================
        if SAVE_GIFS 
            saveGIFs(folder_save);
        end
    end
    
    %% =============== Assign depth non assigned strokes ==================
    
%     assignRemainingStrokes(...
%                         strokes_topology,...
%                         intersections,...
%                         cam_param,...
%                         pairsInterInter);
%                     
% 
%     assignRemainingStrokesHighScore(...
%                         strokes_topology,...
%                         intersections,...
%                         cam_param,...
%                         pairsInterInter);
    
%    assignRemainingStrokesHighScoreBestSolutions(...
%                     strokes_topology,...
%                     intersections,...
%                     cam_param,...
%                     pairsInterInter)
                    
end



function saveGIFs(folder_save)
    global method_str;
    
    filename_gif        = fullfile(folder_save, method_str, 'reconstruction_time.gif');
    filename_gif_depth  = fullfile(folder_save, method_str, 'reconstruction_depth.gif');
    filename_gif_2D     = fullfile(folder_save, method_str, 'reconstruction_time2D.gif');
    %% ============= Save to gif 3D ===================================
    figure(3);
    drawnow;
    frame = getframe(gcf);

    [imind,cm] = rgb2ind(frame.cdata,256);
    if ind_assign == 2 
      imwrite(imind,cm,filename_gif,'gif', 'Loopcount', inf); 
    else 
      imwrite(imind,cm,filename_gif,'gif','WriteMode','append'); 
    end 

    %% ============= Save to gif 3D ===================================
    figure(5);
    drawnow;
    frame = getframe(gcf);

    [imind,cm] = rgb2ind(frame.cdata,256);
    if ind_assign == 2 
      imwrite(imind,cm,filename_gif_depth,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename_gif_depth,'gif','WriteMode','append'); 
    end 

    %% ============= Save to gif 2D ===================================
    figure(4);
    drawnow;
    frame = getframe(gcf);

    [imind,cm] = rgb2ind(frame.cdata,256);
    if ind_assign == 2 
        imwrite(imind,cm,filename_gif_2D,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename_gif_2D,'gif','WriteMode','append'); 
    end 
end