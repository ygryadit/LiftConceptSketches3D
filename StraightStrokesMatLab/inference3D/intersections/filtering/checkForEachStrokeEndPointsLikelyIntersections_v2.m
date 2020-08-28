function mask_intrscts_likely_new = checkForEachStrokeEndPointsLikelyIntersections_v2(intersections, strokes_topology, img)
    %Check for the 25% of the end points if there are likely intersections

    global debug_local;
    debug_local = true;
    global dist_endpoint;
    global SHOW_FIGS_PREPROCESS;
    global USE_CURVES

    
%     intersections_to_pair_new = intersections;
    mask_intrscts_likely_new = false(length(intersections.likely),1);
    
%     for li = 1:length(strokes_topology)
%         %Skip non-straight: 
%         if ~USE_CURVES & strokes_topology(li).primitive_type ~= 0
%             continue;
%         end
%         
%         %Skip non-isolated: 
%         if isempty(strokes_topology(li).indcs_intrsctns)
%               strokes_topology(li).primitive_type = -2;
%               continue;
%         end
%         
%             
%         stroke = strokes_topology(li);
%         stroke.ind = li;
%         
%         %Plot stroke:
%         if debug_local
%             plot_stroke(strokes_topology, li, img);            
%         end
%           
%           
%         % Select intersections with straight stroke only:
%         if ~USE_CURVES
%             mask_intersections_lines = cat(1,strokes_topology(intersections.strokes_indices(stroke.indcs_intrsctns,1)).primitive_type) == 0 & ...
%                                        cat(1,strokes_topology(intersections.strokes_indices(stroke.indcs_intrsctns,2)).primitive_type) == 0;  
% 
%             stroke.indcs_intrsctns = stroke.indcs_intrsctns(mask_intersections_lines);
%             stroke.indcs_intrsctng_strks = stroke.indcs_intrsctng_strks(mask_intersections_lines);
%         end
%         
%         %1. All intersections with previous strokes:
%         mask_prev_strks = stroke.indcs_intrsctng_strks < stroke.ind;
%         indcs_intrsctns_prev_strks = stroke.indcs_intrsctns(mask_prev_strks);
%         
%         if strokes_topology(li).primitive_type == 1
%             %curve
%             intrsctns_seg_nums = stroke.inter_seg_nums(mask_prev_strks);
%         end
%         
%         inter_coord2D_prev_strks = intersections.coordinates2D(indcs_intrsctns_prev_strks, :);
% 
%         %2. Non-collinear intersections with previous strokes:
% %         mask_non_collinear = ~intersections.collinear(indcs_intrsctns_prev_strks);
% %         
% %         ind_non_collinear = indcs_intrsctns_prev_strks(mask_non_collinear);
% %         inter_coord2D_non_collinear = intersections.coordinates2D(ind_non_collinear,:);
% %         
% %         if strokes_topology(li).primitive_type == 1
% %             %curve
% %             intrsctns_seg_nums_non_collinear = intrsctns_seg_nums(mask_non_collinear);
% %         end
%         
% 
%         %3. Likely (based on valency) intersections with previous straight strokes:
% %         ind_likely = indcs_intrsctns_prev_strks(intersections.likely(indcs_intrsctns_prev_strks));
% %         inter_coord2D_likely = intersections.coordinates2D(ind_likely, :);
% %         if debug_local
% %             plot(inter_coord2D_likely(:,1),inter_coord2D_likely(:,2),'c*')
% %         end
%         %% Find distance from each of the end point:
%         length_stroke = stroke.length2D;
%           
%           
% %           stroke_begin = [stroke.points2D(1).x stroke.points2D(1).y];
% %           stroke_end   = [stroke.points2D(end).x stroke.points2D(end).y];
% 
% %         if isempty(inter_coord2D_non_collinear)            
% %             inter_coord2D_non_collinear = inter_coord2D_prev_strks;
% %             if isempty(inter_coord2D_non_collinear)
% %                 continue;
% %             end
% %             if strokes_topology(li).primitive_type == 1
% %                intrsctns_seg_nums_non_collinear = intrsctns_seg_nums;
% %             end                  
% %         end
% 
%         % Find the first and last intersecion along the stroke, assuming they are sorted: 
%         
% %         stroke_begin = inter_coord2D_non_collinear(1,:);
% %         stroke_end   = inter_coord2D_non_collinear(end,:);
% 
%         ind_likely = stroke.indcs_intrsctns(intersections.likely(stroke.indcs_intrsctns));
%         inter_coord2D_likely = intersections.coordinates2D(ind_likely, :);
% 
%         inter_coord2D = intersections.coordinates2D(stroke.indcs_intrsctns, :);
% %         
%         if (size(inter_coord2D_likely,1) < 2)
%            inter_coord2D_likely = inter_coord2D; 
%         end
%     
%         if isempty(inter_coord2D_prev_strks) |  ...
%            (size(inter_coord2D,1) < 2)
%            continue; 
%         end
%         
%      
% 
%         stroke_begin = inter_coord2D_likely(1,:);
%         stroke_end   = inter_coord2D_likely(end,:);
%         
%         if strokes_topology(li).primitive_type == 0
%             distance_begin = sqrt(sum((inter_coord2D_prev_strks - repmat(stroke_begin, size(inter_coord2D_prev_strks,1), 1)).^2,2))./length_stroke;
%             distance_end   = sqrt(sum((inter_coord2D_prev_strks - repmat(stroke_end, size(inter_coord2D_prev_strks,1), 1)).^2,2))./length_stroke;
%         else
%             [distance_begin, distance_end] = computeIntersectionsCurveEndPointsDistance(intrsctns_seg_nums, stroke);
%         end
%         
%         %Intersection near stroke begin:
%         inds_intrscts_likely_new =  ...
%           selectIntersectionsNearExtrem_v2(indcs_intrsctns_prev_strks, ...
%                                            distance_begin,...
%                                            intersections);
%         mask_intrscts_likely_new(inds_intrscts_likely_new) = true;
% 
%         %Intersection near stroke end:
%         inds_intrscts_likely_new =  ...
%           selectIntersectionsNearExtrem_v2(indcs_intrsctns_prev_strks, ...
%                                            distance_end,...
%                                            intersections);
%         mask_intrscts_likely_new(inds_intrscts_likely_new) = true;
%           
%         if li == 31
%             disp('');
%         end
%     end
    
    %% Additional likely intersections for curves:
    if USE_CURVES
       % Find all tangential intersections 
       mask_intrscts_likely_new = mask_intrscts_likely_new | intersections.tangent;
       
       % Find all non-collinear curve-curve
       va = [strokes_topology(intersections.strokes_indices(:,1)).primitive_type] == 1;
       vb = [strokes_topology(intersections.strokes_indices(:,2)).primitive_type] == 1;
       curve_curve = va & vb;
       mask_curve_curve_non_collinear = ~intersections.collinear & curve_curve';
       mask_intrscts_likely_new = mask_intrscts_likely_new | mask_curve_curve_non_collinear;
    end
    
    
    %% Additional intersections:
    %     Finally, if there are strokes that have nothing within 25% of the
    %     end check if there can be added some intersections, if not remove
    %     the stroke from computations.
    
    for li = 1:length(strokes_topology)
        %Skip non-straight: 
        if (strokes_topology(li).primitive_type ~= 0) | ...
            (USE_CURVES & (strokes_topology(li).primitive_type ~= 1) )
            continue;
        end
        
        %Skip non-isolated: 
        if isempty(strokes_topology(li).indcs_intrsctns)
              strokes_topology(li).primitive_type = -2;
              continue;
        end
        
        stroke = strokes_topology(li);
        stroke.ind = li;
        
        %Plot stroke:
        if debug_local
            close(figure(1)); 
            figure(1); 
            imshow(img);
            plot_stroke(strokes_topology, li, img);            
            
        end
                    
        % Select intersections with straight stroke only:
        if (strokes_topology(li).primitive_type == 0)
            mask_intersections_lines = cat(1,strokes_topology(intersections.strokes_indices(stroke.indcs_intrsctns,1)).primitive_type) == 0 & ...
                                       cat(1,strokes_topology(intersections.strokes_indices(stroke.indcs_intrsctns,2)).primitive_type) == 0;  

            stroke.indcs_intrsctns = stroke.indcs_intrsctns(mask_intersections_lines);
            stroke.indcs_intrsctng_strks = stroke.indcs_intrsctng_strks(mask_intersections_lines);
        end
        
        % Find all likely intersections: 
        mask_likely = mask_intrscts_likely_new(stroke.indcs_intrsctns);
        indcs_intrsctns_likely = stroke.indcs_intrsctns(mask_likely);
        
        
        %Plot stroke:
        if debug_local
            close(figure(1)); 
            figure(1); 
            imshow(img);
            plot_stroke(strokes_topology, li, img);            
            plot(intersections.coordinates2D(indcs_intrsctns_likely,1), ...
                 intersections.coordinates2D(indcs_intrsctns_likely,2),'b*');
            strokes_inds = intersections.strokes_indices(indcs_intrsctns_likely,:);
            strokes_inds = setdiff(strokes_inds(:), li);
            for i = strokes_inds'
                 plot(cat(1, strokes_topology(i).points2D(:).x),cat(1, strokes_topology(i).points2D(:).y), 'g');       
            end
        end
        
        
        if isempty(stroke.indcs_intrsctns)
              strokes_topology(li).primitive_type = -2;
              continue;
        end
        
        % Find the first and last intersecion along the stroke, assuming they are sorted: 
        stroke_begin = intersections.coordinates2D(stroke.indcs_intrsctns(1),:);
        stroke_end   = intersections.coordinates2D(stroke.indcs_intrsctns(end),:);
           
        inter_coord2D_likely = intersections.coordinates2D(indcs_intrsctns_likely,:);        
        length_stroke = stroke.length2D;
        
        distance_begin_likely = sqrt(sum((inter_coord2D_likely - repmat(stroke_begin, size(inter_coord2D_likely,1), 1)).^2,2))./length_stroke;
        distance_end_likely   = sqrt(sum((inter_coord2D_likely - repmat(stroke_end, size(inter_coord2D_likely,1), 1)).^2,2))./length_stroke;
        
       
        % Begin:        
        mask_intrscts_likely_new = ...
                assignAdditionalIntersectionsNearExtrem(intersections,... 
                                                        stroke,...
                                                        mask_intrscts_likely_new,...
                                                        stroke_begin,...
                                                        indcs_intrsctns_likely,...
                                                        distance_begin_likely,...
                                                        dist_endpoint,...
                                                        length_stroke,...
                                                        SHOW_FIGS_PREPROCESS,...
                                                        debug_local);
        
        % End:   
        mask_intrscts_likely_new = ...
                assignAdditionalIntersectionsNearExtrem(intersections,... 
                                                        stroke,...
                                                        mask_intrscts_likely_new,...
                                                        stroke_end,...
                                                        indcs_intrsctns_likely,...
                                                        distance_end_likely,...
                                                        dist_endpoint,...
                                                        length_stroke,...
                                                        SHOW_FIGS_PREPROCESS,...
                                                        debug_local);
                                                    
                                                    
%         intersection_near_end_likely = indcs_intrsctns_likely(distance_end_likely < dist_endpoint);
%         
%         if isempty(intersection_near_end_likely )
%             mask_non_collinear = ~intersections.collinear(stroke.indcs_intrsctns);            
%             indcs_intrsctns_nc = stroke.indcs_intrsctns(mask_non_collinear);
%             
%             inter_coord2D = intersections.coordinates2D(indcs_intrsctns_nc,:);
%             distance_end_all   = sqrt(sum((inter_coord2D - repmat(stroke_end, size(inter_coord2D,1), 1)).^2,2))./length_stroke;
%             
%             intersection_near_end_all = indcs_intrsctns_nc(distance_end_all < dist_endpoint);
%             mask_intrscts_likely_new(intersection_near_end_all) = true;
%             
%             if SHOW_FIGS_PREPROCESS || debug_local
%                 plot(intersections.coordinates2D(intersection_near_end_all, 1),...
%                     intersections.coordinates2D(intersection_near_end_all, 2), '*m');    
%             end
%         end

        
    end

end


function mask_intrscts_likely_new = ...
                assignAdditionalIntersectionsNearExtrem(intersections,... 
                                                        stroke,...
                                                        mask_intrscts_likely_new,...
                                                        stroke_begin,...
                                                        indcs_intrsctns_likely,...
                                                        distance_begin_likely,...
                                                        dist_endpoint,...
                                                        length_stroke,...
                                                        SHOW_FIGS_PREPROCESS,...
                                                        debug_local)
    % Likely intersections:
    intersection_near_begin_likely = indcs_intrsctns_likely(distance_begin_likely < dist_endpoint);
        
    if isempty(intersection_near_begin_likely )
        
        mask_non_collinear = ~intersections.collinear(stroke.indcs_intrsctns);

        indcs_intrsctns_nc = stroke.indcs_intrsctns(mask_non_collinear);
        inter_coord2D = intersections.coordinates2D(indcs_intrsctns_nc,:);
        
        mask_likely = intersections.likely(indcs_intrsctns_nc);

        distance_begin_all = sqrt(sum((inter_coord2D - repmat(stroke_begin, size(inter_coord2D,1), 1)).^2,2))./length_stroke;

        intersection_near_begin_likely = indcs_intrsctns_nc(distance_begin_all(mask_likely) < dist_endpoint);

        if ~isempty(intersection_near_begin_likely)
            mask_intrscts_likely_new(intersection_near_begin_likely) = true;
        else
            intersection_near_begin_all = indcs_intrsctns_nc(distance_begin_all < dist_endpoint);
            mask_intrscts_likely_new(intersection_near_begin_all) = true;
            if SHOW_FIGS_PREPROCESS || debug_local
                plot(intersections.coordinates2D(intersection_near_begin_all, 1),...
                    intersections.coordinates2D(intersection_near_begin_all, 2), '*m');    
            end
        
            if isempty(intersection_near_begin_all)
                indcs_intrsctns_c = stroke.indcs_intrsctns(~mask_non_collinear);
                inter_coord2D = intersections.coordinates2D(indcs_intrsctns_c,:);
                distance_begin = sqrt(sum((inter_coord2D - repmat(stroke_begin, size(inter_coord2D,1), 1)).^2,2))./length_stroke;
                intersection_near_begin = indcs_intrsctns_c(distance_begin < dist_endpoint);
                mask_intrscts_likely_new(intersection_near_begin) = true;
                if SHOW_FIGS_PREPROCESS || debug_local
                    plot(intersections.coordinates2D(intersection_near_begin, 1),...
                        intersections.coordinates2D(intersection_near_begin, 2), 'om');    
                end
            end
        end

        
    end

end

function plot_stroke(strokes_topology, li, img)
    figure(14);
    hold off;
    imshow(img);
    hold on;    
    for i = 1:(li-1)
        plot(cat(1, strokes_topology(i).points2D(:).x),cat(1, strokes_topology(i).points2D(:).y), 'b:');       
    end
    plot(cat(1, strokes_topology(li).points2D(:).x),cat(1, strokes_topology(li).points2D(:).y), 'r');       
end

