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
% -------------------------------------------------------------------------
% Output:
% -------------------------------------------------------------------------
% 
%   intersections = 
% 
%       struct with fields:
% 
%       coordinates2D: [num_intersections×2 double] coordinates of the
%                       intersections.
%        line_indices: [num_intersections×2 double] indices only among the
%                       strokes that were considered to be nearly straight
%     p_dist_str_segs: [num_intersections×1 double]
%                       probablity based on distances from the intersection
%                       to each of strokes (unused?)
%     strokes_indices: [num_intersections×2 double]
%                       indices among al the strokes in the
%                       strokes_topology data structure
%           collinear: [num_intersections×1 logical]
%                       is true if the two stokes are either angularily
%                       similar or were marked as lines going towards the
%                       same vanishing point
% 
% -------------------------------------------------------------------------
%   ind_strokes_lines1:
%   ind_strokes_lines2:
%       the indices of all lines that intersect even if they intersect
%       outside the sketching area.
% 
% -------------------------------------------------------------------------
%   ind_do_not_intersect:
%       the indices of those strokes in  [ind_strokes_lines1;
%       ind_strokes_lines2] that ieither intersect outside sketchign are or
%       the intersection point is too far from the stoke end points.


function [  intersections, ...
            ind_strokes_lines1, ...
            ind_strokes_lines2, ...
            ind_do_not_intersect] = ...
                computeIntersectionsBetweenLinePrimitives(...
                        strokes_topology,...
                        h,...
                        w,...
                        img)
global SHOW_FIGS_PREPROCESS;

%% Compute intersections between lines:
ind_lines = find(cat(1, strokes_topology(:).primitive_type) == 0);     

lines = cat(1, strokes_topology(ind_lines).primitive_geom);

lines_groups = cat(1, strokes_topology(ind_lines).line_group);

accuracy_radiuses = cat(1, strokes_topology(ind_lines).accuracy_radius);

[points2D, ind_lines1, ind_lines2] = computeLinesIntersections(lines);


%% Remove points outside image boundaries:
ind = find(points2D(:,1) > 0 & points2D(:,1) < h & ...
        points2D(:,2) > 0 & points2D(:,2) < w);                           

% Fill in intersections structure:
intersections.coordinates2D = points2D(ind,:);
intersections.line_indices  = [ind_lines1(ind) ind_lines2(ind)];

%% Keep only the intersections in some neighborhood around strokes end points.
% Rearrange lines indices, so that they are ordered how they were drawn:
ind_rearrange = find(intersections.line_indices(:,1) > intersections.line_indices(:,2));

temp = intersections.line_indices(ind_rearrange,1);
intersections.line_indices(ind_rearrange,1) = intersections.line_indices(ind_rearrange,2);
intersections.line_indices(ind_rearrange,2) = temp;


lines1_coordinates  = lines(intersections.line_indices(:,1), :);
lines2_coordinates  = lines(intersections.line_indices(:,2), :);

% Distances from the strokes segments to the intersection point:
[distances1, ~] = findLineSegmentPointDistance(lines1_coordinates, intersections.coordinates2D);
[distances2, ~] = findLineSegmentPointDistance(lines2_coordinates, intersections.coordinates2D);    


% for i = 1:length(intersections.coordinates2D)
%     
%     figure(3);    
%     hold off;
%     imshow(img);
%     hold on;    
%     
%     plot(lines(intersections.line_indices(i,1), [1,2]), lines(intersections.line_indices(i,1), [3,4]));
%     plot(lines(intersections.line_indices(i,2), [1,2]), lines(intersections.line_indices(i,2), [3,4]));
%     plot(intersections.coordinates2D(i,1), intersections.coordinates2D(i,2), '*');  
%     legend(sprintf('d1 %d', distances1(i)), sprintf('d2 %d', distances2(i)), 'int');
% end




strokes_thresholds_1  = accuracy_radiuses(intersections.line_indices(:,1), :);
strokes_thresholds_2  = accuracy_radiuses(intersections.line_indices(:,2), :);

% thresholds_intersections = min(accuracy_radiuses(intersections.line_indices), [], 2);



% Probablility based on the length of the strokes:
% p_dist_str_segs = probProximity(distances1, lines1_lengths*thr_length_fraction).*...
%                   probProximity(distances2, lines2_lengths*thr_length_fraction);
%            

% Probablility based on the length of the speed:
% sqrt(2*log(0.5) -- strokes_thresholds_2 is mapped to 0.5 probability

% p_dist_str_segs = probProximity(distances1, strokes_thresholds_1 / sqrt(2*abs(log(0.5))) ).*...
%                   probProximity(distances2, strokes_thresholds_2 / sqrt(2*abs(log(0.5))) );
              
p_dist_str_segs = probProximity(distances1, strokes_thresholds_2 / sqrt(2*abs(log(0.5))) ).*...
                  probProximity(distances2, strokes_thresholds_2 / sqrt(2*abs(log(0.5))) );

ind_intersect = find(p_dist_str_segs > 0.1);


intersections.coordinates2D    = intersections.coordinates2D(ind_intersect, :);
intersections.line_indices  = intersections.line_indices(ind_intersect, :);
intersections.p_dist_str_segs  = p_dist_str_segs(ind_intersect);
intersections.strokes_indices = ind_lines(intersections.line_indices); 

num_pairs = 1:size(ind_lines1,1);
ind_intersect = ind(ind_intersect);
ind_do_not_intersect = [setxor(num_pairs, ind_intersect)];

%% Find which lines are likely to be oversketched:
ind_colinear = findCoincidingLines(intersections.line_indices, lines, img);

lines_groups_s1 = lines_groups(intersections.line_indices(:,1));
lines_groups_s2 = lines_groups(intersections.line_indices(:,2));
% 
% length(ind_colinear)
ind_colinear = [ind_colinear;  find((lines_groups_s1 == lines_groups_s2) & (lines_groups_s2 ~= 4))];
% ind_colinear = [ind_colinear;];
% length(ind_colinear)
% 
% colchck=[lines(intersections.line_indices(:,1),1) lines(intersections.line_indices(:,1),3) ...
%          lines(intersections.line_indices(:,2),1) lines(intersections.line_indices(:,2),3) ...
%          lines(intersections.line_indices(:,2),2) lines(intersections.line_indices(:,2),4)];
% colchck=colchck(:,1).*(colchck(:,4)-colchck(:,6))+colchck(:,3).*(colchck(:,6)-colchck(:,2))+...
%     colchck(:,5).*(colchck(:,2)-colchck(:,4));
% 
% ind_non_colinear=find(abs(colchck)>50);
% intersections.direction_types = false(size(intersections.coordinates2D,1),1);
% intersections.direction_types(ind_non_colinear) = true;

intersections.collinear = false(size(intersections.coordinates2D,1),1);
intersections.collinear(ind_colinear) = true;

ind_strokes_lines1 = ind_lines(ind_lines1);
ind_strokes_lines2 = ind_lines(ind_lines2);

%% Exclude intersection near vanishing points:
% Exclude intersections that are the intersection between strokes of the
% same line group but are not collinear:

mask_same_group = ( (lines_groups_s1 == lines_groups_s2) & ...
                    (lines_groups_s2 ~= 4) & ...
                    (lines_groups_s2 ~= 5) );
                
impossible3D_intersections = (~intersections.collinear & mask_same_group);

intersections.coordinates2D = intersections.coordinates2D(~impossible3D_intersections, :);
intersections.line_indices = intersections.line_indices(~impossible3D_intersections, :);
intersections.p_dist_str_segs = intersections.p_dist_str_segs(~impossible3D_intersections, :);
intersections.strokes_indices = intersections.strokes_indices(~impossible3D_intersections, :);
intersections.collinear = intersections.collinear(~impossible3D_intersections, :);

if SHOW_FIGS_PREPROCESS
    close all;
    hf = figure; imshow(img);
    hold on;
    plot(intersections.coordinates2D(~intersections.collinear,1),...
    intersections.coordinates2D(~intersections.collinear,2),...
    '*');
    plot(intersections.coordinates2D(intersections.collinear,1),...
    intersections.coordinates2D(intersections.collinear,2),...
    'o');
    % legend()

    for i =1:length(lines)
        plot(lines(i,[1,2]),lines(i,[3,4]));
    end

    set(hf, 'Name', 'Approximate intersections');
end

end


