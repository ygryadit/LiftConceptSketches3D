function [intersections,...
    strokes_topology,...
    inds_pairs_collinear_straight_strokes] = computeAllIntersections(strokes_topology, h, w, img)

global SHOW_FIGS_PREPROCESS;
global fid; 
global ACCOUNT_INTER_CURVS;
%% Compute intersections between line primitives:
tic;
[intersections_line_primitives, ...
  ind_strokes_lines1, ...
  ind_strokes_lines2,...
  ind_do_not_intersect] = computeIntersectionsBetweenLinePrimitives(strokes_topology, h, w, img);

elapsedTime = toc;
fprintf(fid, 'Intial estimate of intersections %.3f\n', elapsedTime);

%% Recompute the intersection points carefully between line segements:
tic
[intersections,strokes_topology] = ...
      recomputeAccurateIntersectionsBetweenLineStrokes(...
            strokes_topology,...
            intersections_line_primitives,...
            img);

elapsedTime = toc;
fprintf(fid, 'Refined estimate of intersections %.3f\n', elapsedTime);

 %% Points of near parallel closely located strokes:
tic
ind_strokes_lines1 = ind_strokes_lines1(ind_do_not_intersect);
ind_strokes_lines2 = ind_strokes_lines2(ind_do_not_intersect);

points_attraction = computeAttractionPointsNearParallelStrokes(...
    strokes_topology, ...
    ind_strokes_lines1, ...
    ind_strokes_lines2, ...
    img);

elapsedTime = toc;
fprintf(fid, 'Points of near parallel closely located strokes %.3f\n', elapsedTime);

%% Merge together attraction points and intersection points:
num_attraction_points = size(points_attraction.coordinates2D,1);

intersections.coordinates2D(end+1:end+num_attraction_points,:)      = points_attraction.coordinates2D;
intersections.strokes_indices(end+1:end+num_attraction_points,:)    = points_attraction.strokes_indices;
intersections.collinear(end+1:end+num_attraction_points,:)          = points_attraction.collinear;
intersections.p_dist_str_segs(end+1:end+num_attraction_points,:)    = points_attraction.p_dist_str_segs;
intersections.seg_nums(end+1:end+num_attraction_points,:)           = points_attraction.seg_nums;

clear('points_attraction');


%% Get list of all pairs of collinear lines:
inds_pairs_collinear_straight_strokes = getPairsCollinearStraightStrokes(intersections);


%% Intersections with curves and between curves:
if ACCOUNT_INTER_CURVS
    intersections_c = intersectionsWithCurves(strokes_topology, img);
    num_ic= size(intersections_c.coordinates2D,1);

    intersections.tangent  = false(size(intersections.coordinates2D,1),1);    
    intersections.coordinates2D(end+1:end+num_ic,:)      = intersections_c.coordinates2D;
    intersections.strokes_indices(end+1:end+num_ic,:)    = intersections_c.strokes_indices;
    intersections.collinear(end+1:end+num_ic,:)          = intersections_c.collinear;
    intersections.tangent(end+1:end+num_ic,:)            = intersections_c.tangent;
    intersections.p_dist_str_segs(end+1:end+num_ic,:)    = intersections_c.p_dist_str_segs;
    intersections.seg_nums(end+1:end+num_ic,:)           = intersections_c.seg_nums;
    clear('intersections_c');
end

%% Visualisations:
if SHOW_FIGS_PREPROCESS
    global folder_save_imgs;
    hf = figure; imshow(img);
    hold on;
    plot(intersections.coordinates2D(~intersections.collinear,1),...
        intersections.coordinates2D(~intersections.collinear,2),...
        'm*');
    plot(intersections.coordinates2D(logical(intersections.collinear),1),...
         intersections.coordinates2D(logical(intersections.collinear),2),...
         'y*');

    set(hf, 'Name', 'Exact intersections');
    legend('non collinear', 'collinear');
    saveas(hf, fullfile(folder_save_imgs, 'exact_intersections.png'));
end
    
end