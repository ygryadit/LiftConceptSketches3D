function intersections = mergeIntersecionsSameStrokes(intersections, img, strokes_topology)

%     intersections.likely = intersections.likely';
    
    ind_swap = find(intersections.strokes_indices(:,1) > intersections.strokes_indices(:,2));
    for i = 1:length(ind_swap)
       temp =  intersections.strokes_indices(ind_swap(i),1);
       intersections.strokes_indices(ind_swap(i),1) = intersections.strokes_indices(ind_swap(i),2);
       intersections.strokes_indices(ind_swap(i),2) = temp;
    end
    
    
    [~,ind_unique,~] = unique([intersections.strokes_indices], 'rows');
    for i = 1:length(ind_unique)
        
        if sum([strokes_topology(intersections.strokes_indices(ind_unique(i),:)).primitive_type] == 1)
            continue;
        end
        
        mask = ismember(intersections.strokes_indices,...
                        intersections.strokes_indices(ind_unique(i),:),...
                        'rows');
        inds_merge = find(mask);
       
       if length(inds_merge) >1 
           inds_pairs = nchoosek(inds_merge, 2);
           distances = sqrt(sum((intersections.coordinates2D(inds_pairs(:,1),:) -...
                        intersections.coordinates2D(inds_pairs(:,2),:)).^2,2));
       else
           distances = 0;
       end
       
       
       
%        figure(10);
%        hold off;
%        imshow(img);
%        hold on;       
% %        
%        plot(cat(1,strokes_topology(intersections.strokes_indices(ind_unique(i),1)).points2D(:).x),...
%             cat(1,strokes_topology(intersections.strokes_indices(ind_unique(i),1)).points2D(:).y));
%        plot(cat(1,strokes_topology(intersections.strokes_indices(ind_unique(i),2)).points2D(:).x),...
%             cat(1,strokes_topology(intersections.strokes_indices(ind_unique(i),2)).points2D(:).y));
%        plot(intersections.coordinates2D(inds_merge,1), intersections.coordinates2D(inds_merge,2), '*');
       
       if intersections.collinear(ind_unique(i))
           intersections.coordinates2D(ind_unique(i),:) = ...
                mean(intersections.coordinates2D(inds_merge,:),1);
%            plot(intersections.coordinates2D(ind_unique(i),1), intersections.coordinates2D(ind_unique(i),2), 'o');
%            a =1;
       else
           if length(inds_merge) > 1  
%                try
%                 lines = cat(1, strokes_topology(intersections.strokes_indices(ind_unique(i),:)).primitive_geom);
%                 intersections.coordinates2D(ind_unique(i),:) = computeLinesIntersections(lines);
%                catch e
%                    rethrow(e);
%                end
                intersections.coordinates2D(ind_unique(i),:) = ...
                    mean(intersections.coordinates2D(inds_merge,:),1);
%                  plot(intersections.coordinates2D(ind_unique(i),1), intersections.coordinates2D(ind_unique(i),2), 'o');
%                 a =1;
           end
       end
        
%        figure(11);
%        hold off;
%        imshow(img);
%        hold on;       
%        
%        plot(strokes_topology(intersections.strokes_indices(ind_unique(i),1)).primitive_geom(1:2),...
%             strokes_topology(intersections.strokes_indices(ind_unique(i),1)).primitive_geom(3:4));
%        plot(strokes_topology(intersections.strokes_indices(ind_unique(i),2)).primitive_geom(1:2),...
%             strokes_topology(intersections.strokes_indices(ind_unique(i),2)).primitive_geom(3:4));
%        plot(intersections.coordinates2D(ind_unique(i),1), intersections.coordinates2D(ind_unique(i),2), '*'); 
        
       intersections.accuracy_radius(ind_unique(i),:) = max(distances) + max(intersections.accuracy_radius(inds_merge,:));
    end

    intersections.collinear = intersections.collinear(ind_unique);
    intersections.coordinates2D = intersections.coordinates2D(ind_unique,:);
    intersections.p_dist_str_segs = intersections.p_dist_str_segs(ind_unique,:);
    intersections.seg_nums = intersections.seg_nums(ind_unique,:);
    intersections.strokes_indices = intersections.strokes_indices(ind_unique,:);
    intersections.accuracy_radius = intersections.accuracy_radius(ind_unique,:);
%     intersections.strokes_indices_original = intersections.strokes_indices_original(ind_unique,:);
    
%     intersections.likely = intersections.likely(ind_unique,:);
end