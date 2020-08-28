% -------------------------------------------------------------------------
% Input:
% -------------------------------------------------------------------------
%   intersections:
% ----------------------------
%       struct with fields:
% 
%           coordinates2D: [num_intersections×2 double] 2D coordinates of
%                          each of the 2D intersection points.
%         strokes_indices: [num_intersections×2 double]
%                          pairs of indices of the interescting
%                          strokes_topology
%               collinear: [num_intersections×1 double]
%                          mask of the intersections that are the
%                          intersections of near collinear strokes.
%         p_dist_str_segs: [num_intersections×1 double]
%                           probablility that the intersection is close to
%                           both end points.
%                seg_nums: [num_intersections×2 double]
%                           polyline position in terms of segments in both
%                           strokes.
%                 tangent: [num_intersections×1 logical]
%             
% ----------------------------
%   strokes_topology:
% ----------------------------
%       num_strokes×1 struct array with fields:
% 
%         points2D
%         primitive_type
%         primitive_geom
%         mean_pressure
%         length2DFull
%         length3D
%         length2DPrimitive
%         line_group
%         speed
%         accuracy_radius
%         poly2d_extended
% ----------------------------
% -------------------------------------------------------------------------
% Output:
% -------------------------------------------------------------------------
% ----------------------------
%   strokes_topology:
% ----------------------------
%   All previous fields plus new ones:
%       inds_intrsctns_eval
%       inds_intrsctns_eval_actv
%       inds_intrsctns_eval_mltpl_cnddts
%
%       inds_intrsctng_strks_eval
%       inds_intrsctng_strks_eval_actv
%       inds_intrsctng_strks_eval_mltpl_cnddts
%             primitive_geom_3D
%                   Rough approximation, e.g. in case of line it is the
%                   first and last points of the line stroke.
%             points3D
%             score
%             score_alignment
%             score_coverage
%             confidence
%             direction_vec
%             depth_assigned
%             num_candidate_lines
%             length2D:
%               length between the first and last intersection in the
%               strokes
 
function [strokes_topology] = fillInStrokesData(strokes_topology, intersections, img, INCLUDE_CURVES)
global USE_CURVES;

for li = 1:length(strokes_topology)
   
    % Intersections:
    
    ind1 = find(intersections.strokes_indices(:,1) == li);
    ind2 = find(intersections.strokes_indices(:,2) == li);
    
    if ~INCLUDE_CURVES
       if strokes_topology(li).primitive_type == 0
          stk_inds1 =  intersections.strokes_indices(ind1,2);
          stk_inds2 =  intersections.strokes_indices(ind2,1);
          ind1 = ind1(cat(1,strokes_topology(stk_inds1).primitive_type) == 0);
          ind2 = ind2(cat(1,strokes_topology(stk_inds2).primitive_type) == 0);
       end        
    end
    
    strokes_topology(li).length2D = [];
    strokes_topology(li).indcs_intrsctng_strks = ...
                                [intersections.strokes_indices(ind1,2); ...
                                 intersections.strokes_indices(ind2,1)];
     
    strokes_topology(li).indcs_intrsctns = [ind1; ind2];

    %% Sort the intersections along the line:
    if strokes_topology(li).primitive_type == 0
        inter_coord2D = intersections.coordinates2D(strokes_topology(li).indcs_intrsctns, :);
        stroke_begin = [strokes_topology(li).points2D(1).x strokes_topology(li).points2D(1).y];
        stroke_end   = [strokes_topology(li).points2D(end).x strokes_topology(li).points2D(end).y];

        [~, ...
          ind_sorted] = getSorted2DIntersectionCoordinatesAlongLineDir(inter_coord2D, stroke_begin, stroke_end);

        strokes_topology(li).indcs_intrsctns          = strokes_topology(li).indcs_intrsctns(ind_sorted);
        strokes_topology(li).indcs_intrsctng_strks    = strokes_topology(li).indcs_intrsctng_strks(ind_sorted);
    else
        %for curve sort indices of segments:
        inter_seg_nums = intersections.seg_nums(strokes_topology(li).indcs_intrsctns, :);
        inter_strks_inds = intersections.strokes_indices(strokes_topology(li).indcs_intrsctns, :);
        mask = ismember(inter_strks_inds(:), li);
        inter_seg_nums = inter_seg_nums(:);
        inter_seg_nums = inter_seg_nums(mask);
        [inter_seg_nums,ind_sorted] = sort(inter_seg_nums);
        strokes_topology(li).indcs_intrsctns          = strokes_topology(li).indcs_intrsctns(ind_sorted);
        strokes_topology(li).indcs_intrsctng_strks    = strokes_topology(li).indcs_intrsctng_strks(ind_sorted);
        strokes_topology(li).inter_seg_nums = inter_seg_nums;
    end
%     %% Find distance from each of the end point:
%     strokes_topology(li).distance_begin = sqrt(sum((inter_coord2D - repmat(stroke_begin, size(inter_coord2D,1), 1)).^2,2));
%     strokes_topology(li).distance_end   = sqrt(sum((inter_coord2D - repmat(stroke_end, size(inter_coord2D,1), 1)).^2,2));
%     
    
    % Intersections hypothesis:
    strokes_topology(li).inds_intrsctns_eval = strokes_topology(li).indcs_intrsctns(strokes_topology(li).indcs_intrsctng_strks < li);   
    strokes_topology(li).inds_intrsctng_strks_eval = strokes_topology(li).indcs_intrsctng_strks(strokes_topology(li).indcs_intrsctng_strks < li);
    %strokes_topology(li).mask_intrsctns_prvs_strks  = strokes_topology(li).indcs_intrsctng_strks < li;   
    
    
   
%     strokes_topology(li).coordinates2D = strokes_coordinates_2D(li,:);
    % Intialise 3D data:
    strokes_topology(li).primitive_geom_3D = NaN*ones(length(strokes_topology(li).primitive_geom)/2, 3);
    strokes_topology(li).points3D = NaN*ones(size(strokes_topology(li).points2D,1),3);
    
    % Scores:
    strokes_topology(li).score = 0;
    strokes_topology(li).score_alignment = 0;
    strokes_topology(li).score_coverage = 0;
%     strokes_topology(li).score_intersections = 0;
    
    % Confidence:
    strokes_topology(li).confidence = 0;
    strokes_topology(li).direction_vec = NaN*ones(1, 3);    
    strokes_topology(li).depth_assigned = false;
     
    strokes_topology(li).num_candidate_lines = 0;
    
    
    %Length between last and first intersection points:
%     try

%     strokes_topology(li).length2D = computeStrokeLengthBetweenIntersections(strokes_topology, li, intersections);
    
    if strokes_topology(li).primitive_type == 0
%         disp('----');
%         disp(strokes_topology(li).length2DPrimitive );
        strokes_topology(li).length2D = computeStrokeLengthBetweenIntersectionsApprox(strokes_topology, li, intersections);
%         disp(strokes_topology(li).length2D );
    else
        if USE_CURVES
            if ~isempty(strokes_topology(li).indcs_intrsctns)
                ind_inter1 = 1;
                ind_inter2 = length(strokes_topology(li).indcs_intrsctns);
                try
                    strokes_topology(li).length2D = computeStrokeLengthBetweenIntersections(strokes_topology(li), ind_inter1, ind_inter2);
                catch e
                    rethrow(e);
                end
            else
                strokes_topology(li).length2D = 0;
                strokes_topology(li).primitive_type = -2;
            end
        else
          strokes_topology(li).length2D = NaN;  
        end
        
%         strokes_topology(li).length2D = NaN;
    end
    
    strokes_topology(li).points3D_clean = [];
%     catch
%         dips(li);
%     end
    
    
%     lengthStroke([cat(1, strokes_topology(i).points2D(:).x) cat(1, strokes_topology(i).points2D(:).y)]);
    
    
%     strokes_topology(li).length2D = sqrt(sum((strokes_topology(li).coordinates2D([1,3]) - strokes_topology(li).coordinates2D([2,4])).^2));
    
%     strokes_topology(li).threshold_upper_distance = NaN;

    %     strokes_topology(li).depth_assigned = false;
    
    
%     close(figure(16));
%     figure(16); imshow(img); hold on;
%     plot(strokes_coordinates_2D(li, [1,2]), strokes_coordinates_2D(li, [3,4]));
%    
%     for j = 1:length(strokes_topology(li).indices_intersecting_lines)
%         lj = strokes_topology(li).indices_intersecting_lines(j);
%         plot(strokes_coordinates_2D(lj, [1,2]), strokes_coordinates_2D(lj, [3,4]));
%         plot(intersections.coordinates2D(strokes_topology(li).indcs_intrsctns(j), 1), ...
%             intersections.coordinates2D(strokes_topology(li).indcs_intrsctns(j), 2),...
%             '*');
%     end

    
    
%       plot(intersections.coordinates2D(strokes_topology(li).indcs_intrsctns, 1), ...
%             intersections.coordinates2D(strokes_topology(li).indcs_intrsctns, 2),...
%             '*');
end

end



% function length2D = computeStrokeLengthBetweenIntersections(strokes_topology, li, intersections)
%     if strokes_topology(li).primitive_type ~= 0 || length(strokes_topology(li).indcs_intrsctns) < 2
%         length2D = strokes_topology(li).length2DFull;
%         return;
%     end
%     
%  
%     
%     
%     mask_stroke_in_intersecting_pairs = (intersections.strokes_indices(strokes_topology(li).indcs_intrsctns,:) == li);
%     
%     
%     seg_nums = intersections.seg_nums(strokes_topology(li).indcs_intrsctns,:)';
%     
%     seg_nums = seg_nums(mask_stroke_in_intersecting_pairs');
%     
%     min_seg_num = min(seg_nums);
%     max_seg_num = max(seg_nums);
%     
%     
%     min_seg_num_l = max(floor(min_seg_num),1);
%     min_seg_num_t = max(ceil(min_seg_num),1);
%     
%     min_dir = [strokes_topology(li).points2D(min_seg_num_t).x - strokes_topology(li).points2D(min_seg_num_l).x;
%                 strokes_topology(li).points2D(min_seg_num_t).y - strokes_topology(li).points2D(min_seg_num_l).y];
% %     min_dir = min_dir/norm(min_dir);
%     
%     num_segments =length(strokes_topology(li).points2D(:))-1;
%     
%     max_seg_num_l = max(min(floor(max_seg_num),num_segments),1);
%     max_seg_num_t = min(ceil(max_seg_num),num_segments+1 );
%     try
%         max_dir = [strokes_topology(li).points2D(max_seg_num_t).x - strokes_topology(li).points2D(max_seg_num_l).x;
%                 strokes_topology(li).points2D(max_seg_num_t).y - strokes_topology(li).points2D(max_seg_num_l).y];
%     catch
%         error('');
%     end
% %     max_dir = max_dir/norm(max_dir);
%     
%     first_x = strokes_topology(li).points2D(min_seg_num_l).x + ...
%               min_dir(1) * (min_seg_num - min_seg_num_l);
%           
%     first_y = strokes_topology(li).points2D(min_seg_num_l).y + ...
%               min_dir(2) * (min_seg_num - min_seg_num_l);          
%           
%     last_x = strokes_topology(li).points2D(max_seg_num_l).x + ...
%               max_dir(1) * (max_seg_num - max_seg_num_l);
%           
%     last_y = strokes_topology(li).points2D(max_seg_num_l).y + ...
%               max_dir(2) * (max_seg_num - max_seg_num_l);          
%           
%    
%     length2D = lengthStroke([ [first_x; cat(1, strokes_topology(li).points2D(min_seg_num_t:max_seg_num_l).x); last_x]...        
%                               [first_y; cat(1, strokes_topology(li).points2D(min_seg_num_t:max_seg_num_l).y); last_y] ]);
%     
%     
% end