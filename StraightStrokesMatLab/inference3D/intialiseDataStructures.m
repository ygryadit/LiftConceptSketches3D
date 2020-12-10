% Output:
%   failed is true when the method can not find a vanishign point:
function [strokes_topology, ...
           intersections, ...
           cam_param,...
           pairsInterInter,...
           ind_first_stroke,...
           failed] = ...
    intialiseDataStructures()

    failed = false;
    global folder_save;
    global SHOW_FIGS_PREPROCESS;
    global fid;
    global ORTHOGRAPHIC;
    
    %% Indtialisize:
     intersections = [];
     cam_param = [];
     pairsInterInter = [];
     ind_first_stroke =[];
    
    %% Read sketch and image:  
    tic;
       
    [sketch, img,accuracy_radiuses] = loadSketch();
    global sketch_height;
  
    
    elapsedTime = toc;
    fprintf(fid, 'Step1: Preprocessing and reading sketch: %.3f seconds\n', elapsedTime);
    
    %% Strokes topology, separate curves vs lines:
    tic;
    
    strokes_topology = findAllStrokeLines(sketch.strokes);
    
    %% Save separation curves strokes:
    save_as_svg_line_strokes(strokes_topology, folder_save, 'lines.svg');
    save_as_svg_curve_strokes(strokes_topology, folder_save, 'curves.svg');
    save_as_svg_sketch(strokes_topology, folder_save, 'sketch.svg');
    
    
    elapsedTime = toc;
    fprintf(fid, 'Step 2: Separation to straight and curved strokes: %.3f seconds\n', elapsedTime);
   
%     thr_mark = max(5,max(cat(1, strokes_topology(:).length2DFull))*0.05);
%     mask_mark = find(cat(1, strokes_topology(:).length2DFull) < thr_mark);
%     
%     for i = mask_mark'  
%         strokes_topology(i).primitive_type = 2;
%     end
%     %% Clip strokes to the drawing area:
% 
%     strokes_topology = trimStraightStrokesCoordinatesDrawingField(strokes_topology);

    %% Estimate VP:
    fprintf(fid, 'Step 3: Vanishing points ...\n');
    tic
    
    [probabilitiesVP, VP,...
     vps_selected,... % additional vanishing points for sets of paralel lines
     inds_axis,...    % axis that is orthogonal to the line
     inds_active_lines ... % indices of lines that belond to a certain cluster
     ] = getVanishingPoints(strokes_topology);
    
    if size(VP,1) < 3
        ORTHOGRAPHIC = true;
        failed = true;
        return;
    end
    
    elapsed_time = toc;
    fprintf(fid, '\t  Elapsed time = %.3f seconds\n', elapsed_time);
    
    %% Rearrange stroke coordinates:
    [~, vp_indices] = max(probabilitiesVP,[],2);                            
   
%     strokes_topology = rearrangeLineCoordinates(strokes_topology, VP, vp_indices);


    
    %% Estimate camera parameters:
 
    fprintf(fid, 'Step 4: Compute camera projection matrix ...');
   
    tic;
    
    [lines_group, VP, vp_new_ind] = ...
                estimateInitialCameraParameters( ...
                    VP,...
                    probabilitiesVP,...
                    strokes_topology);
    
                
    new_vp_inds(1) = find(vp_new_ind == 1);
    new_vp_inds(2) = find(vp_new_ind == 2);
    new_vp_inds(3) = find(vp_new_ind == 3);
    
    inds_axis = updateIndOrthAx(inds_axis, new_vp_inds);
      
 
    if SHOW_FIGS_PREPROCESS
        global folderVP;
        mask_lines = cat(1, strokes_topology(:).primitive_type) == 0;
        lines_coordinates2D = cat(1,strokes_topology(mask_lines).primitive_geom);
        visualizeVP(SHOW_FIGS_PREPROCESS,folderVP,img,VP',probabilitiesVP,lines_coordinates2D);
    end
       clear('probabilitiesVP');
    toc;
    elapsed_time = toc;
    fprintf(fid, '\t  Elapsed time = %.3f seconds\n', elapsed_time);
    
    
    %% Add line group:
    strokes_topology = assignLineDirection(strokes_topology, ...
                                           lines_group,...
                                           vps_selected,... % additional vanishing points for sets of paralel lines
                                           inds_axis,...    % axis that is orthogonal to the line
                                           inds_active_lines ... % indices of lines that belond to a certain cluster
                                           );
    
    %% Compute accuracy radiuses based on speed of the stroke:    
    strokes_topology = assignAccuracyRadius(strokes_topology, accuracy_radiuses);
    
    if SHOW_FIGS_PREPROCESS
        plotAccuracyRadiuses(img, strokes_topology);
    end

    
    %% Intersections:
    
    sketch = rmfield(sketch,'strokes');
    
    
    fprintf(fid, 'Step 5: Find a set of valid intersections in 3D:\n');
   
    tic
    [intersections, ...
     pairsInterInter, ...
     strokes_topology] = ...
            findAllIntersectionPointsOfAllStrokeSegements(strokes_topology, ...
                                                          sketch.canvas.height, ...
                                                          sketch.canvas.width, ...
                                                          img);
    elapsed_time = toc;
    fprintf(fid, '\t Step 5: completed in %.3f seconds: \n', elapsed_time);
    
    %% Select the first stroke:
    ind_first_stroke = findFirstSuitableStroke(strokes_topology, intersections);
    
    %% Re-estimate camera parameters so that the camera origin is in the end point of the selected stroke:
    % The strokes were merged so the stroke origin could have changed, plus
    % the first stroke is selected differently:
    point2DConterpartOrigin = ...
        strokes_topology(ind_first_stroke).primitive_geom(1,[1,3]);
    
    [cam_param, ~, vp_new_ind] = estimateCameraParameters(VP, ...
                                                          sketch_height, ...
                                                          point2DConterpartOrigin,...
                                                          img);
                                                      
    cam_param.VP = vp_new_ind;
    cam_param.vp_coord = VP(vp_new_ind,:);

    new_vp_inds(1) = find(vp_new_ind == 1);
    new_vp_inds(2) = find(vp_new_ind == 2);
    new_vp_inds(3) = find(vp_new_ind == 3);
    new_vp_inds(4) = 4; %line with unknown direction
    new_vp_inds(5) = 5; %line in one of the main planes
    new_vp_inds(6) = 6; %curve
    
    
    vals = num2cell(arrayfun(@(x) new_vp_inds(x.line_group), strokes_topology));
    [strokes_topology(:).line_group] = vals{:};
    

    
    vals = num2cell(arrayfun(@(x) updateIndOrthAx(x.ind_orth_ax, new_vp_inds), strokes_topology));
    [strokes_topology(:).ind_orth_ax] = vals{:};
    

%     p_ = probabilitiesVP;
%     p_(:, 1:3) = p_(:, vp_new_ind);
%     
%     % global ZUP;
%     % ZUP = vp_(3,2) < 0.5*(vp_(1,2)+vp_(2,2));
%     
%     [~, lines_group] = max(p_,[],2);    
%     
%     
    
    % Add line group:
%     inds_axis = vp_new_ind(inds_axis);
%     strokes_topology = assignLineDirection(strokes_topology, ...
%                                            lines_group,...
%                                            vps_selected,... % additional vanishing points for sets of paralel lines
%                                            inds_axis,...    % axis that is orthogonal to the line
%                                            inds_active_lines ... % indices of lines that belond to a certain cluster
%                                            );
    
    
    save_as_svg_straigt_strokes_vp_colorcoded(strokes_topology);
    %% Mark short lines as marks:

    strokes_topology = selectMarks(strokes_topology);
    
%     save_as_svg_sketch_radiuses(strokes_topology, folder_save, 'sketch_strokes_radiuses.svg');
    
end


    
function ind_orth_ax = updateIndOrthAx(ind_orth_ax, new_vp_inds)
    if ~isempty(ind_orth_ax) & ~isnan(ind_orth_ax)
        ind_orth_ax = new_vp_inds(ind_orth_ax);
    else
        ind_orth_ax = NaN;
    end
end