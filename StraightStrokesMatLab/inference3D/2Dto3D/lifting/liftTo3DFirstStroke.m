function [stroke, intersections] = liftTo3DFirstStroke(stroke, intersections, cam_param, img)
    stroke.primitive_geom_3D(1,:) = [0,0,0];
    p_line_prior_3D_2 = estimateOneCoordinate(stroke.line_group, ...
                                             stroke.primitive_geom(1,[2,4]), ...
                                             [0,0,0], cam_param.P, img);
    
    % Line second point and polyline points:
    points_2D = [stroke.primitive_geom(1,[2,4]); ...
                 [cat(1, stroke.points2D(:).x), cat(1, stroke.points2D(:).y);...
                  cat(1,intersections(stroke.indcs_intrsctns).coordinates2D) ]];
    
    try
    [points3D, lambdas] = opt3Dpos2DProj(points_2D,...
                                         cam_param, ...
                                         stroke.primitive_geom_3D(1,:), ...
                                         p_line_prior_3D_2);
    catch e
        rethrow(e)
    end
    stroke.primitive_geom_3D(2,:)= points3D(1,:);
    num_points = length(stroke.points2D);
    stroke.points3D = points3D(2:(num_points+1),:);
    stroke.points3D_clean = stroke.points3D;
    
    vals_assign = num2cell(points3D((num_points+2):end,:),2);
    [intersections(stroke.indcs_intrsctns).coordinates3D] = vals_assign{:};
    
%     intersections.coordinates3D(stroke.indcs_intrsctns, :) =  points3D((num_points+2):end,:);
                                         
    stroke.score           = costPairwise(0,0,1,false);    
    stroke.score_alignment = 1.0;
    stroke.confidence      = 1.0;
    stroke.direction_vec   = stroke.primitive_geom_3D(2,:)-stroke.primitive_geom_3D(1,:);
    stroke.length3D = norm(stroke.direction_vec);    
    stroke.direction_vec   = stroke.direction_vec./norm(stroke.direction_vec);
    stroke.depth_assigned = true;
     
%     plotEstimates(stroke,intersections, p_line_prior_3D_2);
%     stroke.threshold_upper_distance = 0.1*sqrt(sum((stroke.coordinates3D(4:6) - stroke.coordinates3D(1:3)).^2));
end



function plotEstimates(stroke,intersections, p_line_prior_3D_2)
    close(figure(1));
    figure(1);
    hold on;
        
    % Plot estimate:
    plot3( [stroke.primitive_geom_3D(:,1)],...
            [stroke.primitive_geom_3D(:,2)],...
            [stroke.primitive_geom_3D(:,3)], 'g*', 'LineWidth', 2);    
        
    % Plot polyline:
    plot3( [stroke.points3D(:,1)],...
            [stroke.points3D(:,2)],...
            [stroke.points3D(:,3)], ':', 'LineWidth', 2);    
     
    c3D = cat(1, intersections(stroke.indcs_intrsctns).coordinates3D); 
    plot3( c3D(:,1),...
            c3D(:,2),...
            c3D(:,3),...
           'ro');       
     
       axis equal;  
end