function intersections_new = checkForEachStrokeEndPointsLikelyIntersections(intersections, strokes_topology, img)
    %Check for the 25% of the end points if there are likely intersections

    global debug_local;
    debug_local = false;
    
    intersections_new = intersections;
    
    for li = 1:length(strokes_topology)
          
          if strokes_topology(li).primitive_type ~= 0
            continue;
          end
          
          if debug_local
              figure(14);
              hold off;
              imshow(img);
              hold on;

              plot(cat(1, strokes_topology(li).points2D(:).x),cat(1, strokes_topology(li).points2D(:).y));                        
          end
          
          
          if isempty(strokes_topology(li).indcs_intrsctns)
              strokes_topology(li).primitive_type = -2;
              continue;
          end
          
          stroke = strokes_topology(li);
          stroke.ind = li;
          
          % Select intersections with straight stroke only:
          mask_intersections_lines = cat(1,strokes_topology(intersections.strokes_indices(stroke.indcs_intrsctns,1)).primitive_type) == 0 & ...
                                     cat(1,strokes_topology(intersections.strokes_indices(stroke.indcs_intrsctns,2)).primitive_type) == 0;  
                              
          stroke.indcs_intrsctns = stroke.indcs_intrsctns(mask_intersections_lines);
          stroke.indcs_intrsctng_strks = stroke.indcs_intrsctng_strks(mask_intersections_lines);
          
          %1.All intersections coordiantes with straihgt strokes:
          inter_coord2D = intersections.coordinates2D(stroke.indcs_intrsctns, :);
          
          %2.Non-collinear intersections with straihgt strokes:
          ind_non_collinear = stroke.indcs_intrsctns(~intersections.collinear(stroke.indcs_intrsctns));
          inter_coord2D_non_collinear = intersections.coordinates2D(ind_non_collinear,:);
          
    
          
          if debug_local
            %3.Likely (based on valency) intersections with straight strokes:
            ind_likely = stroke.indcs_intrsctns(intersections.likely(stroke.indcs_intrsctns));          
            strks_ind_likely = stroke.indcs_intrsctng_strks(intersections.likely(stroke.indcs_intrsctns));          
            
            inter_coord2D_likely = intersections.coordinates2D(ind_likely, :);
            plot(inter_coord2D_likely(:,1),inter_coord2D_likely(:,2), '*r');            
            for mi = reshape(strks_ind_likely, 1, [])
                plot(cat(1, strokes_topology(mi).points2D(:).x),cat(1, strokes_topology(mi).points2D(:).y), 'b:');                 
            end
          end
          
          %% Find distance from each of the end point:
          length_stroke = stroke.length2D; %computeLengthStroke([cat(1,stroke.points2D(:).x), cat(1,stroke.points2D(:).y)]);
          
          
%           stroke_begin = [stroke.points2D(1).x stroke.points2D(1).y];
%           stroke_end   = [stroke.points2D(end).x stroke.points2D(end).y];

          if isempty(inter_coord2D_non_collinear)
             ind_all = stroke.indcs_intrsctns;
             inter_coord2D_non_collinear = intersections.coordinates2D(ind_all,:);
             if isempty(inter_coord2D_non_collinear)
                continue;
             end
          end

          % Find the first and last intersecion along the stroke, assuming they are sorted: 
          stroke_begin = inter_coord2D_non_collinear(1,:);
          stroke_end   = inter_coord2D_non_collinear(end,:);
            

          distance_begin = sqrt(sum((inter_coord2D - repmat(stroke_begin, size(inter_coord2D,1), 1)).^2,2))./length_stroke;
          distance_end   = sqrt(sum((inter_coord2D - repmat(stroke_end, size(inter_coord2D,1), 1)).^2,2))./length_stroke;
          
 
          [intersections_to_pair_new_likely] = selectIntersectionsNearExtrem(stroke, distance_begin, intersections);
          intersections_new.likely = intersections_new.likely | intersections_to_pair_new_likely;
          [intersections_to_pair_new_likely] = selectIntersectionsNearExtrem(stroke, distance_end, intersections);
          intersections_new.likely = intersections_new.likely | intersections_to_pair_new_likely;
          
          % If there are no intersection with preceeding strokes add them:
          mask_likely = intersections_new.likely(stroke.indcs_intrsctns);
          indcs_intrsctng_strks_likely = stroke.indcs_intrsctng_strks(mask_likely);
          
          if isempty(indcs_intrsctng_strks_likely < stroke.ind)
             mask_previous_strokes = stroke.indcs_intrsctng_strks < stroke.ind;
             inds_intersections_activate = stroke.indcs_intrsctns(mask_previous_strokes);
             intersections_new.likely(inds_intersections_activate) = true;
             inter_coord2D_new = intersections_new.coordinates2D(inds_intersections_activate,:);
             if debug_local
                plot(inter_coord2D_new(:,1),inter_coord2D_new(:,2), '*m');
             end
          end
          
          
          
%           intersection_near_stroke_begin = strokes_topology(li).indcs_intrsctns(distance_begin < dist_endpoint);
%           intersection_near_stroke_begin = intersection_near_stroke_begin(~intersections.collinear(intersection_near_stroke_begin));
%           
%           
%           if isempty(intersection_near_stroke_begin) %first intersection is very far from the stroke begin point:
%               intersection_near_stroke_begin = strokes_topology(li).indcs_intrsctns(1);
%           end
%           
%           intersection_near_stroke_end = strokes_topology(li).indcs_intrsctns(distance_end < dist_endpoint);
%           intersection_near_stroke_end = intersection_near_stroke_end(~intersections.collinear(intersection_near_stroke_end));
%           
%           if isempty(intersection_near_stroke_end)
%               intersection_near_stroke_end = strokes_topology(li).indcs_intrsctns(end);
%           end
%           
% %           plot(intersections.coordinates2D( strokes_topology(li).indcs_intrsctns, 1),...
% %                 intersections.coordinates2D( strokes_topology(li).indcs_intrsctns, 2), 'ob'); 
%           
%          
% 
%           if ~sum(intersections.likely(intersection_near_stroke_begin))% & sum(intersections.likely(intersection_near_stroke_end))
%              %No intersection near the end point is identified as likely:
%              intersections.likely(intersection_near_stroke_begin) = true;
%              if SHOW_FIGS_PREPROCESS
%                  plot(intersections.coordinates2D(intersection_near_stroke_begin, 1),...
%                     intersections.coordinates2D(intersection_near_stroke_begin, 2), 'og'); 
%              end
%           else
%             if debug_local
%                 intersection_near_stroke_begin_likely = intersection_near_stroke_begin(intersections.likely(intersection_near_stroke_begin));
%                 plot(intersections.coordinates2D(intersection_near_stroke_begin_likely, 1),...
%                     intersections.coordinates2D(intersection_near_stroke_begin_likely, 2), 'or'); 
%             end
%           end
%           
%           
%           if ~sum(intersections.likely(intersection_near_stroke_end))% & sum(intersections.likely(intersection_near_stroke_begin))
%              %No intersection near the end point is identified as likely:
%              intersections.likely(intersection_near_stroke_end) = true;
%              if SHOW_FIGS_PREPROCESS
%                   plot(intersections.coordinates2D(intersection_near_stroke_end, 1),...
%                     intersections.coordinates2D(intersection_near_stroke_end, 2), 'og'); 
%              end
%          else
%             if debug_local
%                 intersection_near_stroke_end_likely = intersection_near_stroke_end(intersections.likely(intersection_near_stroke_end));
%                 plot(intersections.coordinates2D(intersection_near_stroke_end_likely, 1),...
%                     intersections.coordinates2D(intersection_near_stroke_end_likely, 2), 'or'); 
%             end
%           end
     end

end


