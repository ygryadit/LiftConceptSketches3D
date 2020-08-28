function [lines_group, VP, vp_new_ind] = ...
                estimateInitialCameraParameters( ...
                    VP,...
                    probabilitiesVP,...
                    strokes_topology)
    global sketch_height;  
    global filepath_sketch_img;     
    global ORTHOGRAPHIC;
        
    img = readSketchImg(filepath_sketch_img);
    
   
    
    if ~ORTHOGRAPHIC
        %% Find first suitable stroke:            
        ind_first_stroke = ...
            findFirstVisibleStrokeTowardsVP(...
                        probabilitiesVP,...
                        strokes_topology);

        point2DConterpartOrigin = strokes_topology(ind_first_stroke).primitive_geom(1,[1,3]);
    
        [cam_param, VP, vp_new_ind] = ...
                estimateCameraParameters(   VP, ...
                                            sketch_height, ...
                                            point2DConterpartOrigin,...
                                            img);
                                        
        cam_param.VP = vp_new_ind;
        cam_param.vp_coord = VP(vp_new_ind,:);
        p_ = probabilitiesVP;
        p_(:, 1:3) = p_(:, vp_new_ind);

        [~, lines_group] = max(p_,[],2);    
    else
        cam_param = NaN;
        lines_group = linesGroupsOrthographicCamera(strokes_topology);        
    end

    
    
    
end

