function sketch_strokes = resampleAllStrokes(sketch_strokes, DGP_THR, img)
    % Resampling constant:
    epsilon = DGP_THR;

    for i = 1:length(sketch_strokes)
%         figure(2);
% % 
%         imshow(img);
%         hold on;
%         plot(curve(i).coordinates2D(:,1), curve(i).coordinates2D(:,2),'ro-');
%         disp(size(curve(i).coordinates2D,1));
        sketch_strokes_points = [cat(1, sketch_strokes(i).points.x) cat(1, sketch_strokes(i).points.y)]';
        if isempty(sketch_strokes_points)
            continue;
        end
        [~,indices] = DouglasPeucker(sketch_strokes_points,1:size(sketch_strokes_points,2),epsilon);
        
        p_vals_original = cat(1,sketch_strokes(i).points(:).p);
        
        sketch_strokes(i).points = sketch_strokes(i).points(indices);
        sketch_strokes_points = sketch_strokes_points';
%         disp(cat(1,sketch_strokes(i).points(:).p));

        sketch_strokes_points_e = sketch_strokes_points(2:end,:);
        sketch_strokes_points_b = sketch_strokes_points(1:end-1,:);
        sq_dist = (sketch_strokes_points_e - sketch_strokes_points_b).^2;
        dist_segs = sqrt(sum( sq_dist,2 ));
        mean_p = sum(p_vals_original(2:end).*dist_segs ,1)./...
                                            sum(dist_segs,1);
        for j = 1:length(indices)
            
            
%             disp(sketch_strokes(i).points(j).p);
            
            sketch_strokes(i).points(j).p = mean_p;
%             disp(sketch_strokes(i).points(j).p);
        end
        
        %Smoothen
        
%          for j = 2:(length(indices)-1)
%              sketch_strokes(i).points(j).p  = mean(p_vals_original(j-1:j+1));
%          end
%          sketch_strokes(i).points(1).p = sketch_strokes(i).points(2).p;
%          sketch_strokes(i).points(end).p = sketch_strokes(i).points(end-1).p;
%          disp(cat(1,sketch_strokes(i).points(:).p));
         
%         plot(sketch_strokes_points(1,:), sketch_strokes_points(2,:), 'r');
%         plot(sketch_strokes_points(1,indices), sketch_strokes_points(2,indices), 'b');


    end

end