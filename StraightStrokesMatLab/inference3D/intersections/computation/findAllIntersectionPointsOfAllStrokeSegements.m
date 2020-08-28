% -------------------------------------------------------------------------
% Input:
% -------------------------------------------------------------------------
%       strokes_topology: 
% 
%           num_strokes×1 struct array with fields:
% 
%           points2D
%           primitive_type
%           primitive_geom
%           mean_pressure
%           length2DFull
%           length3D
%           length2DPrimitive
%           line_group
%           speed
%           accuracy_radius
% 
%       
%       h,w:
%           dimensions of the drawing fields
% 
% -------------------------------------------------------------------------
% Output:
% -------------------------------------------------------------------------
% 
%       intersections -- structure with the followong fields:
% 
%           coordinates2D:
%               2D coordinates of each of the 2D intersection points.
%           strokes_indices:
%               pairs of indices of the interescting strokes
%           p_dist_str_segs:
%               probablility that the intersection is close to both end points
%           direction_types:
%               false if the strokes are likely to be alligned, true if the
%               intersection is likely to be a true intersection.
% 
%       points_alligned -- structure with the followong fields:
% 
%          coordinates2D:
%               2D coordinates of each of the 2D points, where the strokes are
%               closely located and are nearly alligned.
%          strokes_indices:
%               pairs of indices of the strokes
%           p_dist_str_segs:
%               probablility that the intersection is close to the end points
%               of strokes.

function [  intersections_, ...            
            pairsInterInter, ...           
            strokes_topology] = ...
                findAllIntersectionPointsOfAllStrokeSegements(strokes_topology, h,w, img)
    global fid;
    global SHOW_FIGS_PREPROCESS;
    global folder_save_imgs;
    
    global ACCOUNT_INTER_CURVS;
    global USE_CURVES;
%     %% Compute intersections between line primitives:
%     tic;
%     [intersections_line_primitives, ...
%       ind_strokes_lines1, ...
%       ind_strokes_lines2,...
%       ind_do_not_intersect] = computeIntersectionsBetweenLinePrimitives(strokes_topology, h, w, img);
%   
%     elapsedTime = toc;
%     fprintf(fid, 'Intial estimate of intersections %.3f\n', elapsedTime);
%     
%     %% Recompute the intersection points carefully between line segements:
%     tic
%     [intersections,strokes_topology] = ...
%           recomputeAccurateIntersectionsBetweenLineStrokes(...
%                 strokes_topology,...
%                 intersections_line_primitives,...
%                 img);
% 
%     elapsedTime = toc;
%     fprintf(fid, 'Refined estimate of intersections %.3f\n', elapsedTime);
%     
%     if SHOW_FIGS_PREPROCESS
%         hf = figure; imshow(img);
%         hold on;
%         plot(intersections.coordinates2D(~intersections.collinear,1),...
%             intersections.coordinates2D(~intersections.collinear,2),...
%             'm*');
%         plot(intersections.coordinates2D(logical(intersections.collinear),1),...
%              intersections.coordinates2D(logical(intersections.collinear),2),...
%              'y*');
%          
%         set(hf, 'Name', 'Exact intersections');
%         legend('non collinear', 'collinear');
%         saveas(hf, fullfile(folder_save_imgs, 'exact_intersections.png'));
%     end
%     %% Points of near parallel closely located strokes:
%     tic
%     ind_strokes_lines1 = ind_strokes_lines1(ind_do_not_intersect);
%     ind_strokes_lines2 = ind_strokes_lines2(ind_do_not_intersect);
%     
%     points_attraction = computeAttractionPointsNearParallelStrokes(...
%         strokes_topology, ...
%         ind_strokes_lines1, ...
%         ind_strokes_lines2, ...
%         img);
% 
%     elapsedTime = toc;
%     fprintf(fid, 'Points of near parallel closely located strokes %.3f\n', elapsedTime);
%     
%     %% Merge together attraction points and intersection points:
%     num_attraction_points = size(points_attraction.coordinates2D,1);
%     
%     intersections.coordinates2D(end+1:end+num_attraction_points,:)      = points_attraction.coordinates2D;
%     intersections.strokes_indices(end+1:end+num_attraction_points,:)    = points_attraction.strokes_indices;
%     intersections.collinear(end+1:end+num_attraction_points,:)          = points_attraction.collinear;
%     intersections.p_dist_str_segs(end+1:end+num_attraction_points,:)    = points_attraction.p_dist_str_segs;
%     intersections.seg_nums(end+1:end+num_attraction_points,:)           = points_attraction.seg_nums;
%     
%     clear('points_attraction');
%     
%     
%     %% Get list of all pairs of collinear lines:
%     inds_pairs_collinear_straight_strokes = getPairsCollinearStraightStrokes(intersections);
%     
%     
%     %% Intersections with curves and between curves:
%     if ACCOUNT_INTER_CURVS
%         intersections_c = intersectionsWithCurves(strokes_topology, img);
%         num_ic= size(intersections_c.coordinates2D,1);
% 
%         intersections.tangent  = false(size(intersections.coordinates2D,1),1);    
%         intersections.coordinates2D(end+1:end+num_ic,:)      = intersections_c.coordinates2D;
%         intersections.strokes_indices(end+1:end+num_ic,:)    = intersections_c.strokes_indices;
%         intersections.collinear(end+1:end+num_ic,:)          = intersections_c.collinear;
%         intersections.tangent(end+1:end+num_ic,:)            = intersections_c.tangent;
%         intersections.p_dist_str_segs(end+1:end+num_ic,:)    = intersections_c.p_dist_str_segs;
%         intersections.seg_nums(end+1:end+num_ic,:)           = intersections_c.seg_nums;
%         clear('intersections_c');
%     end
   
    
    [intersections,...
   strokes_topology,...
    inds_pairs_collinear_straight_strokes] = computeAllIntersections(strokes_topology, h, w, img);

    inds_pairs_collinear_strokes = getPairsCollinearStraightStrokes(intersections);

    
%     %% Find connected pairs of intersections:
    intersections = ...
        assignAccuracyRadiusLatestStroke(strokes_topology, intersections);
    
    % plotAccuracyRadiusIntersections(intersections, img);
    intersections = reorderIntersectons(intersections);
%     
%     [pairsInterInter] = ...
%                     findConnectedPairsOfIntersections(strokes_topology, intersections, img);


    %% Fill in stroke topology structure non merged strokes:
    INCLUDE_CURVES = true;
    [strokes_topology] = fillInStrokesData(strokes_topology, intersections, img, INCLUDE_CURVES);
    
    object_str = jsonencode(intersections);
%     fid_ = fopen(fullfile(folder_save, ...
%                      sprintf('%s_%s_%s.json', ...
%                      designer,...
%                      object_name,...
%                      'intersections_before_merging')), ...
%                      'w');
%     fwrite(fid_, object_str);
%     fclose(fid_); 
    
    %% Merge lines:
    save_as_svg_primitives(strokes_topology, 'lines_before_merging.svg');
    global folder_save;
    save_as_svg_curve_strokes(strokes_topology, fullfile(folder_save, 'merged_svg'), 'curves_before_merging.svg');
    
    if SHOW_FIGS_PREPROCESS
        figure;
        hold on;
        for li = 1:length(strokes_topology)
            plot(cat(1,strokes_topology(li).points2D.x),...
                 cat(1,strokes_topology(li).points2D.y));   
        end
    end
    
    
    
    global MERGE_LINES;
    if MERGE_LINES
        for i = 1:length(strokes_topology)
            strokes_topology(i).merged_with = [];
        end
        [intersections, strokes_topology] = ...
                groupCollinearStrokes(inds_pairs_collinear_straight_strokes,...
                                   intersections,...
                                   strokes_topology,...
                                   img);
         if USE_CURVES                      
            [~, strokes_topology] = ...
                    groupCollinearCurvedStrokes(inds_pairs_collinear_strokes,...
                                       intersections,...
                                       strokes_topology,...
                                       img);                        
         end
    end 
    
 

     save_as_svg_primitives(strokes_topology, 'lines_after_merging.svg');
%     save_as_svg_clusters(strokes_topology, 'lines_after_merging_clusters.svg');                                         
%      save_as_svg_curve_strokes(strokes_topology_temp, fullfile(folder_save, 'merged_svg'), 'curves_after_merging.svg'); 
     clear('strokes_topology_temp');
%     [intersections,...
%     strokes_topology,...
%     inds_pairs_collinear_straight_strokes] = computeAllIntersections(strokes_topology, h, w, img);

                                              
%     %% Remove intersections with curves:
%     inds_keep = find(( cat(1,strokes_topology_(intersections.strokes_indices(:,1)).primitive_type) == 0 & ...
%                        cat(1,strokes_topology_(intersections.strokes_indices(:,2)).primitive_type) == 0 ));
%     
%     
%     intersections.coordinates2D      = intersections.coordinates2D(inds_keep,:);
%     intersections.strokes_indices    = intersections.strokes_indices(inds_keep,:);
%     intersections.collinear          = intersections.collinear(inds_keep);    
%     intersections.p_dist_str_segs    = intersections.p_dist_str_segs(inds_keep);
%     intersections.seg_nums           = intersections.seg_nums(inds_keep,:);
%     intersections.likely             = intersections.likely(inds_keep);
%     intersections.accuracy_radius    = intersections.accuracy_radius(inds_keep);

    %% Merge intersections the same stroke
     [intersections,...
     strokes_topology,...
     ~] = computeAllIntersections(strokes_topology, h, w, img);
 
     intersections = ...
        assignAccuracyRadiusLatestStroke(strokes_topology, intersections);
    
    % plotAccuracyRadiusIntersections(intersections, img);
    intersections = reorderIntersectons(intersections);
    
    intersections = mergeIntersecionsSameStrokes(intersections, img, strokes_topology);
    
        
     %% Along each line find which intersections are likely to be reliable:
    inds_pairs_collinear_strokes = getPairsCollinearStraightStrokes(intersections);

    INCLUDE_CURVES = true;
     
    
    intersections = ...
        assignAccuracyRadiusLatestStroke(strokes_topology, intersections);
    
    % plotAccuracyRadiusIntersections(intersections, img);
    intersections = reorderIntersectons(intersections);
    
    [strokes_topology] = fillInStrokesData(strokes_topology, intersections, img, INCLUDE_CURVES);
    
    global sketch_height;
    
    save_as_svg_all_strokes_intersections_all_v2(strokes_topology, intersections, folder_save_imgs, sketch_height, 'all_intersections_all_strokes.svg');
    
    [pairsInterInter] = ...
                    findConnectedPairsOfIntersections(strokes_topology, intersections, img);
                    
    intersections = ...
                findAlongLineRelibleIntersections(strokes_topology,...
                                                  intersections,...
                                                  pairsInterInter,... 
                                                  inds_pairs_collinear_strokes,...
                                                  img);
                                              
                                              
    %% Fill in the strokes_topology after the strokes and intersections are merged:
%     INCLUDE_CURVES = false;
     
%     [intersections, pairsInterInter] = ...
%                     findConnectedPairsOfIntersections(strokes_topology, intersections, img);
                
%     [strokes_topology] = fillInStrokesData(strokes_topology, intersections, img, INCLUDE_CURVES);
    [strokes_topology] = computeLength2DLikely(strokes_topology, intersections);

    
  
    if SHOW_FIGS_PREPROCESS
        figure(15);
        imshow(img); hold on;
        m1 = intersections.likely;
        plot(intersections.coordinates2D(logical(m1),1), intersections.coordinates2D(logical(m1),2),'*');
    end

    %% Copy
    num_intersections = size(intersections.coordinates2D,1);
    for i = 1:num_intersections
        intersections_(i).coordinates2D     = intersections.coordinates2D(i,:);
        intersections_(i).strokes_indices   = intersections.strokes_indices(i,:);
%         intersections_(i).seg_nums          = intersections.seg_nums(i,:);
        intersections_(i).p_dist_str_segs   = intersections.p_dist_str_segs(i,:);
        intersections_(i).accuracy_radius   = intersections.accuracy_radius(i,:);
        intersections_(i).collinear         = intersections.collinear(i,:);
        intersections_(i).likely            = intersections.likely(i);
        intersections_(i).coordinates3D     = NaN*ones(1,3);
        intersections_(i).is_active         = NaN;
    end
    
    
    %% Curves
     global sketch_height;
%      save_as_svg_straigt_strokes_intersections_all(strokes_topology, intersections_, folder_save_imgs, sketch_height, 'all_intersections.svg');
%      save_as_svg_straigt_strokes_intersections_likely(strokes_topology, intersections_, folder_save_imgs, sketch_height, 'likely_intersections.svg');
     
     num_intersections = length(intersections_);
     num_likely_intersections = sum([intersections_(:).likely]);
     save(fullfile(folder_save, 'intersections_stat.mat'), 'num_intersections', 'num_likely_intersections');
     
     
end


