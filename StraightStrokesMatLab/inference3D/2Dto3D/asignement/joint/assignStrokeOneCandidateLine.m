function [stroke, intersections] = assignStrokeOneCandidateLine(stroke, intersections, cam_param,candidate_line)
    ic3D = cat(1,intersections(stroke.indcs_intrsctns).coordinates3D);
    if ~isempty(ic3D)
        ind_stroke_undefined_intersections = stroke.indcs_intrsctns(isnan(ic3D(:, 1)));
    else
        ind_stroke_undefined_intersections = [];
    end
%     if sum(abs(candidate_line.coordinates3D_prior(4:6) - candidate_line.coordinates3D_prior(1:3))) < 1e-10
%         %[TODO::] estimate for one ens, should be careful which one:
% %         stroke.primitive_geom(1,[1,3]); stroke.primitive_geom(1,[2,4]);
%         deg_fredom = 1;
%         points_2D = [stroke.primitive_geom(1,[2,4]); ... % line primitive end points
%                     [cat(1, stroke.points2D(:).x), cat(1, stroke.points2D(:).y)];... %polyline points
%                      intersections.coordinates2D(ind_stroke_undefined_intersections, :) ]; %non_assigned_intersectons
%         
%         
%     else
    deg_fredom = 0;
    if isfield(stroke, 'points2DOriginal') & ~isempty(stroke.points2DOriginal)
        num_points = length(stroke.points2DOriginal);
        points_2D       = [[cat(1, stroke.points2DOriginal(:).x), cat(1, stroke.points2DOriginal(:).y)];... %polyline points
                  cat(1, intersections(ind_stroke_undefined_intersections).coordinates2D) ]; %non_assigned_intersectons
        num_points_clean = length(stroke.points2D);      
        points_2D_clean = [[cat(1, stroke.points2D(:).x), cat(1, stroke.points2D(:).y)];... %polyline points
                  cat(1, intersections(ind_stroke_undefined_intersections).coordinates2D) ]; %non_assigned_intersectons
              
    else
        num_points = length(stroke.points2D);
        points_2D = [[cat(1, stroke.points2D(:).x), cat(1, stroke.points2D(:).y)];... %polyline points
                  cat(1, intersections(ind_stroke_undefined_intersections).coordinates2D) ]; %non_assigned_intersectons
    end
    

    [points3D, lambdas] = opt3Dpos2DProj(points_2D, cam_param, ...
                                         candidate_line.coordinates3D_prior(1, 1:3), ...
                                         candidate_line.coordinates3D_prior(1, 4:6));
    if exist('points_2D_clean', 'var')
        [points3D_clean, lambdas] = opt3Dpos2DProj(points_2D_clean, cam_param, ...
                                         candidate_line.coordinates3D_prior(1, 1:3), ...
                                         candidate_line.coordinates3D_prior(1, 4:6));
    else
        points3D_clean = [];
    end
    
    if deg_fredom
        line_direction =  candidate_line.coordinates3D(1:3) - points3D(1,:); %one intersection is fixed
        %propogate from the interseection point:
        stroke.primitive_geom_3D(1,1:3) = find3DCoordinatesLine(  cam_param.P,...
                                                                candidate_line.coordinates3D(1:3),...
                                                                line_direction,...
                                                                stroke.primitive_geom(1,[1,3]));
        %assign optimised value:
        stroke.primitive_geom_3D(2,:) = points3D(1,:);
    else
        line_direction =  candidate_line.dir; %two intersections are fixed
        %propogate from the interseection point 1:
        stroke.primitive_geom_3D(1,1:3) = candidate_line.coordinates3D_prior(1:3);
%         find3DCoordinatesLine(  cam_param.P,...
%                                                                 candidate_line.coordinates3D(1:3),...
%                                                                 line_direction,...
%                                                                 stroke.primitive_geom(1,[1,3]));
        %propogate from the interseection point 2:
        stroke.primitive_geom_3D(2,1:3) = candidate_line.coordinates3D_prior(4:6);
% find3DCoordinatesLine(  cam_param.P,...
%                                                                 candidate_line.coordinates3D(4:6),...
%                                                                 line_direction,...
%                                                                 stroke.primitive_geom(1,[2,4]));
    end
    
    
    stroke.points3D = points3D((1+deg_fredom):(num_points+deg_fredom),:);
    if ~isempty(points3D_clean)
        stroke.points3D_clean = points3D_clean((1+deg_fredom):(num_points_clean+deg_fredom),:);
    end
    assign_vals = num2cell(points3D((num_points+deg_fredom+1):end,:),2);    
    [intersections(ind_stroke_undefined_intersections).coordinates3D] = assign_vals{:};
%     disp(intersections.coordinates3D(ind_stroke_undefined_intersections, :));                                   
    stroke.score           = candidate_line.max_cost;

    try
        ind = find(cat(1,candidate_line.configurations(:).p_full_joint) == stroke.score);
        stroke.score_alignment = candidate_line.configurations(ind).p_directional; %changes
    catch e
        rethrow(e);
    end
    
    stroke.confidence      = 1.0; %change
%     stroke.direction_vec   = stroke.primitive_geom_3D(2,:)-stroke.primitive_geom_3D(1,:);
    stroke.direction_vec   = line_direction./norm(line_direction);
    stroke.depth_assigned = true;
    stroke.length3D = candidate_line.length3D;
    if isfield(stroke, 'candidate_lines')
        stroke.candidate_lines = [];
    end
    stroke.num_candidate_lines =0;
end