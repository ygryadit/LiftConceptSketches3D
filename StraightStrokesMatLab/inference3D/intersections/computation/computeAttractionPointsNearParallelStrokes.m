function points_attraction = computeAttractionPointsNearParallelStrokes(strokes_topology, ind_strokes_lines1, ind_strokes_lines2, img)

lines1 = cat(1, strokes_topology(ind_strokes_lines1).primitive_geom);
lines2 = cat(1, strokes_topology(ind_strokes_lines2).primitive_geom);



% [distances(:,1), intersection_coordinates(:, 1, :)] = ...
%         find2DLineSegmentPointDistance(lines1, lines2(:,[1,3]));
% 
% [distances(:,2), intersection_coordinates(:, 2, :)] = ...
%         find2DLineSegmentPointDistance(lines1, lines2(:,[2,4]));
%     
% [distances(:,3), intersection_coordinates(:, 3, :)] = ...
%         find2DLineSegmentPointDistance(lines2, lines1(:,[1,3]));
%     
% [distances(:,4), intersection_coordinates(:, 4, :)] = ...
%         find2DLineSegmentPointDistance(lines2, lines1(:,[2,4]));

start_points1 = zeros(length(ind_strokes_lines1), 2);
end_points1 = zeros(length(ind_strokes_lines1), 2);
start_points2 = zeros(length(ind_strokes_lines1), 2);
end_points2 = zeros(length(ind_strokes_lines1), 2);
for i =1:length(ind_strokes_lines1)
    start_points1(i,:)   = [(strokes_topology(ind_strokes_lines1(i)).points2D(1).x) (strokes_topology(ind_strokes_lines1(i)).points2D(1).y)];
    end_points1(i,:)     = [(strokes_topology(ind_strokes_lines1(i)).points2D(end).x) (strokes_topology(ind_strokes_lines1(i)).points2D(end).y)];
    start_points2(i,:) = [(strokes_topology(ind_strokes_lines2(i)).points2D(1).x) (strokes_topology(ind_strokes_lines2(i)).points2D(1).y)];
    end_points2(i,:) = [(strokes_topology(ind_strokes_lines2(i)).points2D(end).x) (strokes_topology(ind_strokes_lines2(i)).points2D(end).y)];
end

[distances(:,1), intersection_coordinates(:, 1, :)] = ...
        find2DLineSegmentPointDistance(lines1, start_points2);

[distances(:,2), intersection_coordinates(:, 2, :)] = ...
        find2DLineSegmentPointDistance(lines1, end_points2);
    
[distances(:,3), intersection_coordinates(:, 3, :)] = ...
        find2DLineSegmentPointDistance(lines2, start_points1);
    
[distances(:,4), intersection_coordinates(:, 4, :)] = ...
        find2DLineSegmentPointDistance(lines2, end_points1);


[best_distances, ind] = min(distances,[],2);

num_attraction_points = size(intersection_coordinates,1); 
points_pairs = zeros(num_attraction_points, 4);

%% Pairs of the points
points_pairs(:,1) = intersection_coordinates(...
                        sub2ind(size(intersection_coordinates),...
                                (1:size(intersection_coordinates,1))',...
                                ind,...
                                ones(size(intersection_coordinates,1), 1)));

points_pairs(:,2) = intersection_coordinates(...
                        sub2ind(size(intersection_coordinates),...
                                (1:size(intersection_coordinates,1))',...
                                ind,...
                                2*ones(size(intersection_coordinates,1), 1)));

points_pairs(:,3) = intersection_coordinates(...
                        sub2ind(size(intersection_coordinates),...
                                (1:size(intersection_coordinates,1))',...
                                ind,...
                                3*ones(size(intersection_coordinates,1), 1)));

points_pairs(:,4) = intersection_coordinates(...
                        sub2ind(size(intersection_coordinates),...
                                (1:size(intersection_coordinates,1))',...
                                ind,...
                                4*ones(size(intersection_coordinates,1), 1)));
                            
%% Compute proximaty criteria
ind_strokes = [ind_strokes_lines1, ind_strokes_lines2];
% ind_stroke_threshold = min([ind_strokes_lines1, ind_strokes_lines2], [], 1);
ind_stroke_threshold = max([ind_strokes_lines1, ind_strokes_lines2], [], 1);

accuracy_radiuses = strokes_topology(ind_stroke_threshold).accuracy_radius;

p_Proximity = probProximity(best_distances, accuracy_radiuses / sqrt(2*abs(log(0.5))) );

ind_attraction = find(p_Proximity > 0.5);

points_attraction.coordinates2D     = points_pairs(ind_attraction,:);
points_attraction.strokes_indices   = ind_strokes(ind_attraction,:);
points_attraction.p_dist_str_segs   = p_Proximity(ind_attraction);
if ~isempty(points_attraction.strokes_indices)
    mask_collinear = checkCollinear(strokes_topology, points_attraction.strokes_indices);
else
    mask_collinear =[];
end

points_attraction.collinear = mask_collinear;

%% Plot
global SHOW_FIGS_PREPROCESS

if SHOW_FIGS_PREPROCESS
    hf = figure; 
    imshow(img);
    hold on;
    for i = 1:length(ind_attraction)
    %     imshow(img);
    %     hold on;
        if (points_attraction.collinear(i))        
            plot(points_attraction.coordinates2D(i,[1,3]), points_attraction.coordinates2D(i,[2,4]), '*-');
        else
            plot(points_attraction.coordinates2D(i,[1,3]), points_attraction.coordinates2D(i,[2,4]), 'o-');
        end
        plot(cat(1,strokes_topology(points_attraction.strokes_indices(i,1)).points2D.x),...
             cat(1,strokes_topology(points_attraction.strokes_indices(i,1)).points2D.y));
        plot(cat(1,strokes_topology(points_attraction.strokes_indices(i,2)).points2D.x),...
             cat(1,strokes_topology(points_attraction.strokes_indices(i,2)).points2D.y));
    %     hold off;
    end
    set(hf, 'Name', 'Attraction points');
end
%Copy individual closest point coordinates:


coordinates2D_strokes = points_attraction.coordinates2D;

points_attraction.seg_nums = convertIntersectionCoord2SegNum(points_attraction.coordinates2D, points_attraction.strokes_indices, strokes_topology);


points_attraction.coordinates2D = [];
%Average x coordinates:
points_attraction.coordinates2D(:,1) = 0.5*(coordinates2D_strokes(:,1) + coordinates2D_strokes(:,3));
%Average y coordinates:
points_attraction.coordinates2D(:,2) = 0.5*(coordinates2D_strokes(:,2) + coordinates2D_strokes(:,4));


end

function seg_nums = convertIntersectionCoord2SegNum(coordinates2D, strokes_indices, strokes_topology)
    % Convertes casterian coordinates to a linear coordiante along the polyline:  
    % Input
    %     coordinates2D N x [x1 y1 x2 y2], where N is a number of attraction
    %       points
    %     strokes_indices
    %     strokes_topology

    % Number of attraction points:
    N = size(coordinates2D,1);
    seg_nums = zeros(N,2);
    for i = 1:N
        seg_nums(i,1) = casterCoord2linearSegCord(coordinates2D(i,1:2), strokes_topology(strokes_indices(i,1)).points2D);
        seg_nums(i,2) = casterCoord2linearSegCord(coordinates2D(i,3:4), strokes_topology(strokes_indices(i,2)).points2D);
    end
end

function seg_num = casterCoord2linearSegCord(coordinates2D, polyline)
   polyline = [cat(1,polyline.x) cat(1,polyline.y)];
   
   distances = sum(polyline - coordinates2D,2).^2;
   [~, inds] = sort(distances);
   
   ind1 = inds(1);
   ind2 = inds(2);
   
   seg_num = ind1 + norm(coordinates2D -  polyline(ind1,:))/norm(polyline(ind2,:) -  polyline(ind1,:));
%    
%    disp(polyline)
%    disp(coordinates2D)
%     disp(seg_num)
end

function mask_collinear= checkCollinear(strokes_topology, strokes_indices)
    lines1 = cat(1,strokes_topology(strokes_indices(:,1)).primitive_geom);
    lines2 = cat(1,strokes_topology(strokes_indices(:,2)).primitive_geom);
    
    dir1 =  lines1(:, [2,4]) - lines1(:, [1,3]);
    dir2 =  lines2(:, [2,4]) - lines2(:, [1,3]);    
    
    dir1 = dir1./repmat(sqrt(sum(dir1.^2,2)), 1, 2);
    dir2 = dir2./repmat(sqrt(sum(dir2.^2,2)), 1, 2);

    cos_dirs = dot(dir1, dir2, 2);


    % disp(acos(cos_dirs)/pi*180)
    mask_collinear = abs(cos_dirs) > cos(pi/9);
end