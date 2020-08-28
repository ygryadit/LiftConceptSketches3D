% Recomputes the accurate positions of the intersections betwee line
% strokes. The positions are inaccurate when the strokes are approximated
% with line segments.
% 
% -------------------------------------------------------------------------
% Input:
% -------------------------------------------------------------------------
%   strokes_topology: 
% -----------------------------
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
% -----------------------------
%   intersections:
% -----------------------------
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


function [intersections_new, strokes_topology] = ...
            recomputeAccurateIntersectionsBetweenLineStrokes(...
                strokes_topology,...
                intersections,...
                img)

    intersections_new.coordinates2D     =[];
    intersections_new.strokes_indices   =[];
    intersections_new.collinear   =[];
    intersections_new.p_dist_str_segs   =[];
    intersections_new.seg_nums = [];
    num_intersections = size(intersections.strokes_indices,1);
    
    strokes_topology = extendStrokes(strokes_topology, img);
     
    for i = 1:num_intersections        
%         if intersections.p_dist_str_segs(i) < (1.0 + 1e-5 % intersections is out of line.
%             %Intersection near the end point, just copy:
%             intersections_new.coordinates2D(end+1,:) =  intersections.coordinates2D(i,:);
%             intersections_new.strokes_indices(end+1,:) = intersections.strokes_indices(i,:);
%             intersections_new.collinear(end+1)=intersections.collinear(i,:);
%             intersections_new.p_dist_str_segs(end+1)=intersections.p_dist_str_segs(i,:);
% %             intersections_new.seg_nums(end+1,:) = [i_i, j_i]; 
%             continue;
%         end
        
        % polylines inidces:
        s1 = intersections.strokes_indices(i,1);
        s2 = intersections.strokes_indices(i,2);
        
        % intersection coordinates:
        
        
        [begin_s1_x, begin_s1_y, end_s1_x, end_s1_y] = extendIntersectingPolyline(strokes_topology(s1));
        [begin_s2_x, begin_s2_y, end_s2_x, end_s2_y] = extendIntersectingPolyline(strokes_topology(s2));
        
        
                                          
        [xi,yi,i_i, j_i] = ...
            intersectionsPolyPoly([begin_s1_x; cat(1,strokes_topology(s1).points2D.x); end_s1_x],...
                                  [begin_s1_y; cat(1,strokes_topology(s1).points2D.y); end_s1_y],...
                                  [begin_s2_x; cat(1,strokes_topology(s2).points2D.x); end_s2_x],...
                                  [begin_s2_y; cat(1,strokes_topology(s2).points2D.y); end_s2_y]);

        
        
        num_ss_intersections = length(xi);
        
        intersections_new.coordinates2D(end+1:end+num_ss_intersections,:) = [xi yi];
        
        intersections_new.seg_nums(end+1:end+num_ss_intersections,:) = [i_i-1, j_i-1]; % coordinates of segements where the intersection is.
        
        intersections_new.strokes_indices(end+1:end+num_ss_intersections,:) = ...
            repmat(intersections.strokes_indices(i,:), num_ss_intersections, 1);
        
%         ib = length(intersections_new.collinear)+1;
        
        intersections_new.collinear(end+1:end+num_ss_intersections)=...
            intersections.collinear(i,1);
        
%         ie = length(intersections_new.collinear);
        
        intersections_new.p_dist_str_segs(end+1:end+num_ss_intersections)=...
            intersections.p_dist_str_segs(i,1);
        
%         
%         figure(3); imshow(img); hold on;
%         plot(cat(1,strokes_topology(s1).points2D.x),cat(1,strokes_topology(s1).points2D.y));
%         plot(cat(1,strokes_topology(s2).points2D.x),cat(1,strokes_topology(s2).points2D.y));        
% %         plot(xi,yi,'*');
%         ind = [ib:ie];
%         ind_col = ind( intersections_new.collinear(ind) == 1);
%         ind_ncol = ind( intersections_new.collinear(ind) == 0);
%         plot(intersections_new.coordinates2D(ind_col,1), intersections_new.coordinates2D(ind_col,2),'o');
%         plot(intersections_new.coordinates2D(ind_ncol,1), intersections_new.coordinates2D(ind_ncol,2),'*');
    end

    intersections_new.collinear = reshape(intersections_new.collinear, [], 1);
    intersections_new.p_dist_str_segs = reshape(intersections_new.p_dist_str_segs, [], 1);
end


