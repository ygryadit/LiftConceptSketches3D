function intersections = ...
                findAlongLineRelibleIntersections(strokes_topology,...
                                            intersections,...
                                            pairsInterInter,... 
                                            collinear_strokes_pairs,...
                                            img)


global SHOW_FIGS_PREPROCESS
global USE_CURVES

if SHOW_FIGS_PREPROCESS
    hf = figure(13);
    hold off;
    imshow(img);
    hold on;
end

inds_corners = findIntersectionPlanesCorners(intersections, strokes_topology);
% intersections.likely(inds_corners) = true;

plot(intersections.coordinates2D(inds_corners(:), 1),...
    intersections.coordinates2D(inds_corners(:), 2), '*g'); 

intersections.likely = false(size(intersections.coordinates2D,1),1);

for int_ind = 1:size(intersections.coordinates2D,1)       
    
   %% Exclude collinear intersections from the set of likely intersections 
   if  intersections.collinear(int_ind)
       intersections.likely(int_ind) = false;
       continue;
   end

   %% Check that the intersections are between two lines
   if ~USE_CURVES
       if  strokes_topology(intersections.strokes_indices(int_ind,1)).primitive_type ~= 0 | ...
           strokes_topology(intersections.strokes_indices(int_ind,2)).primitive_type ~= 0
           intersections.likely(int_ind) = false;
           continue;
       end
   end
   
%     strokes_indices = intersections.strokes_indices(int_ind,:);

%          if SHOW_FIGS_PREPROCESS
%               hf = figure(13);
%                 hold off;
%                 imshow(img);
%                 hold on;
%                 plot(cat(1,strokes_topology(strokes_indices(1)).points2D(:).x),...
%                  cat(1,strokes_topology(strokes_indices(1)).points2D(:).y));
%             
%                 plot(cat(1,strokes_topology(strokes_indices(2)).points2D(:).x),...
%                  cat(1,strokes_topology(strokes_indices(2)).points2D(:).y));
%          end

% 
% if SHOW_FIGS_PREPROCESS
%     hf = figure(13);
%     hold off;
%     imshow(img);
%     hold on;
% end
    if intersections.likely(int_ind)
        continue;
    end
    
   [valency, ind_paired_intersections] = forEachStrokeIntersectionEstimateValence(  int_ind, ...
                                                        pairsInterInter,...
                                                        intersections,...                
                                                        collinear_strokes_pairs,...
                                                        strokes_topology);



    if valency >= 3            
        if SHOW_FIGS_PREPROCESS

%                 plot(cat(1,strokes_topology(strokes_indices(1)).points2D(:).x),...
%                  cat(1,strokes_topology(strokes_indices(1)).points2D(:).y));
%             
%                 plot(cat(1,strokes_topology(strokes_indices(2)).points2D(:).x),...
%                  cat(1,strokes_topology(strokes_indices(2)).points2D(:).y));

            plot(intersections.coordinates2D(int_ind, 1),...
                intersections.coordinates2D(int_ind, 2), '*r'); 
            plot(intersections.coordinates2D(ind_paired_intersections, 1),...
                intersections.coordinates2D(ind_paired_intersections, 2), 'or'); 
        end
%             text(intersections.coordinates2D(int_ind, 1),...
%                 intersections.coordinates2D(int_ind, 2), num2str(int_ind))
%             disp('');
        intersections.likely(int_ind) = true;
        intersections.likely(ind_paired_intersections) = true;
    else
        intersections.likely(int_ind) = false;
    end
end

inds_intersections_corners = find(intersections.likely);

%% Keep likely only withing 25% from th eend with previous stroke, include collinear if needed:
% intersections.planar(inds_corners(:)) = true;

intersections = ...
    checkForEachStrokeEndPointsLikelyIntersections(...
        intersections,...
        strokes_topology,...
        img);
    

global folder_save_imgs;
global sketch_height;


save_as_svg_straigt_strokes_intersections_all_v2(strokes_topology, intersections, folder_save_imgs, sketch_height, 'all_intersections_straight.svg');
save_as_svg_curved_strokes_intersections_all_v2(strokes_topology, intersections, folder_save_imgs, sketch_height, 'all_intersections_curved.svg');
save_as_svg_straigt_strokes_intersections_likely_v2(strokes_topology, intersections, inds_intersections_corners, folder_save_imgs, sketch_height, 'likely_intersections.svg');
     

% intersections.likely = ...
%     checkForEachStrokeEndPointsLikelyIntersections_v2(...
%         intersections,...
%         strokes_topology,...
%         img);
    
%% Close planar loops:
% mask = intersections.likely(inds_corners);
% inds_rows = find(sum(mask,2) == 3);
% inds_intrsctns_change_to_likely = inds_corners(inds_rows,:);
% inds_intrsctns_change_to_likely = unique(inds_intrsctns_change_to_likely(mask(inds_rows,:) == 0));
% intersections.likely(inds_intrsctns_change_to_likely) = true;

% for attract_ind = 1:size(points_attraction_to_pair.coordinates2D,1)    
%       valency = forEachStrokeAttractionEstimateValence(attract_ind, ...
%                         intersections,...
%                         pairsInterAttrac,...
%                         points_attraction_to_pair,...
%                         collinear_strokes_pairs,...
%                         img, strokes_topology);
%                     
% %         figure(13);
% %         hold off;
% %         imshow(img);
% %         hold on;
% %         i = intersections.strokes_indices(int_ind,1);
% %         plot(cat(1,strokes_topology(i).points2D(:).x), cat(1,strokes_topology(i).points2D(:).y));
% %         i = intersections.strokes_indices(int_ind,2);
% %         plot(cat(1,strokes_topology(i).points2D(:).x), cat(1,strokes_topology(i).points2D(:).y));
% %         
%         if valency >= 3
% %             plot(cat(1,strokes_topology(i).points2D(:).x), cat(1,strokes_topology(i).points2D(:).y));
% %             plot(points_attraction_to_pair.coordinates2D(attract_ind, [1,3]),...
% %                 points_attraction_to_pair.coordinates2D(attract_ind, [2,4]), '*b');  
% % 
% %             points_attraction_to_pair.likely(attract_ind) = true;
%         else
%             points_attraction_to_pair.likely(attract_ind) = false;
%         end
% end

% for i = 1:length(strokes_topology)
% %     all_relevant_stroeks = find(i == strokes_quaternions(:,1) | ...
% %          i == strokes_quaternions(:,2) | ...
% %          i == strokes_quaternions(:,3) | ...
% %          i == strokes_quaternions(:,4));
%     stroke_intersection_inds =  find((intersections.strokes_indices(:,1) == i) | ...
%                                     (intersections.strokes_indices(:,2) == i) );
%                                 
%     for j = 1:length(stroke_intersection_inds)
%        %For each intersection find all paired intersections
%        int_ind = stroke_intersection_inds(j);
%        
%        
%        valency = forEachStrokeIntersectionEstimateValence(int_ind, ...
%                     pairsInterInter,...
%                     intersections,...
%                     pairsInterAttrac,...
%                     points_attraction_to_pair,...
%                     collinear_strokes_pairs, img);
%        
%              
%         if valency >= 3
%             plot(cat(1,strokes_topology(i).points2D(:).x), cat(1,strokes_topology(i).points2D(:).y));
%             plot(intersections.coordinates2D(int_ind, 1),...
%                 intersections.coordinates2D(int_ind, 2), '*r');  
% %             disp('');
%             intersections.likely(int_ind) = true;
%         else
%             intersections.likely(int_ind) = false;
%         end
%     end
%     
%     
%     
%     stroke_attraction_inds =  find((points_attraction_to_pair.strokes_indices(:,1) == i) | ...
%                                    (points_attraction_to_pair.strokes_indices(:,2) == i) );
%                                
%     for j = 1:length(stroke_attraction_inds)                          
%          attract_ind = stroke_attraction_inds(j);
%          valency = forEachStrokeAttractionEstimateValence(attract_ind, ...
%                         intersections,...
%                         pairsInterAttrac,...
%                         points_attraction_to_pair,...
%                         collinear_strokes_pairs,...
%                         img);
%         if valency >= 3
%             plot(cat(1,strokes_topology(i).points2D(:).x), cat(1,strokes_topology(i).points2D(:).y));
%             plot(points_attraction_to_pair.coordinates2D(attract_ind, [1,3]),...
%                 points_attraction_to_pair.coordinates2D(attract_ind, [2,4]), '*b');  
% 
%             points_attraction_to_pair.likely(attract_ind) = true;
%         else
%             points_attraction_to_pair.likely(attract_ind) = false;
%         end
%         
%     end
% end

if SHOW_FIGS_PREPROCESS
    global folder_save_imgs;
    saveas(hf, fullfile(folder_save_imgs, 'intersections_likely.png'));
end
end



%  selectIntersectionsNearEndPoint(strokes_topology(li), distances_extrem, intersections)




function intersections = ...
            addSecondCollinearIntersection(intersections, ...
                        ind_intrsctn)

  
%        if intersections.collinear(ind_intrsctn)
%             % Find the other intersection with the collinear stroke:
%             ind_inter_add = find(cat(1,intersections.strokes_indices(:,1)) == ...
%                                     intersections.strokes_indices(ind_intrsctn,1),...
%                                  cat(1,intersections.strokes_indices(:,2)) == ...
%                                     intersections.strokes_indices(ind_intrsctn,2));
%                                 
%             intersections.likely(ind_inter_add) = true;                   
%        end
end



function  valency = forEachStrokeAttractionEstimateValence(attract_ind, ...
                        intersections,...
                        pairsInterAttrac,...
                        points_attraction_to_pair,...
                        collinear_strokes_pairs,...
                        img, strokes_topology)

     %% Intersections:
       ind_paired_intersections = pairsInterAttrac((pairsInterAttrac(:,2) == attract_ind),1);
       if ~isempty(ind_paired_intersections) 
           ind_lines_intersect = [points_attraction_to_pair.strokes_indices(attract_ind,:);...
               intersections.strokes_indices([ind_paired_intersections],:)];
       
            
%             figure(13);
%             imshow(img);
%             hold on;
      
           
           valency = estimateValenceLinesGroup(ind_lines_intersect, collinear_strokes_pairs,  intersections.strokes_indices, strokes_topology);
%            
%            if valency >=3
%                 for i = ind_lines_intersect
%                     plot(cat(1,strokes_topology(i).points2D(:).x), cat(1,strokes_topology(i).points2D(:).y));
%                 end
%            end
       else
           valency = 0;
           disp("That has to be checked");
       end
end

function [valency, ind_paired_intersections] = forEachStrokeIntersectionEstimateValence(int_ind, ...
                                                            pairsInterInter,...
                                                            intersections,...
                                                            collinear_strokes_pairs,...
                                                            strokes_topology)
       
       %% Intersections:
       ind_paired_intersections1 = pairsInterInter((pairsInterInter(:,1) == int_ind),2);
       ind_paired_intersections2 = pairsInterInter((pairsInterInter(:,2) == int_ind),1);

       
       
       ind_paired_intersections = [ind_paired_intersections1; ind_paired_intersections2];   

       
       
%         plot(intersections.coordinates2D(ind_paired_intersections, 1),...
%              intersections.coordinates2D(ind_paired_intersections, 2), '*b'); 
       
       ind_lines_intersect = intersections.strokes_indices([int_ind; ind_paired_intersections],:);
       
        
       %% Estimate valence:
       
%        ind_lines_intersect = unique(ind_lines_intersect(:))';
       ind_lines_intersect = unique(ind_lines_intersect, 'rows');
       ind_lines_intersect = orderPairs(ind_lines_intersect);
       
%        maks_swap = find(ind_lines_intersect(:,2) < ind_lines_intersect(:,1));
%        temp = ind_lines_intersect(maks_swap, 2);
%        ind_lines_intersect(maks_swap, 2) = ind_lines_intersect(maks_swap, 1);
%        ind_lines_intersect(maks_swap, 1) = temp; 
       
%        
%         figure(13);
%         imshow(img);
%         hold on;
       ind_lines_intersect = unique(sort(ind_lines_intersect(:)));
       
       
%        global SHOW_FIGS_PREPROCESS
% 
%         if SHOW_FIGS_PREPROCESS
%            for i = 1:length(ind_lines_intersect)
%                plot([strokes_topology(ind_lines_intersect(i)).points2D(:).x], ...
%                     [strokes_topology(ind_lines_intersect(i)).points2D(:).y]);
%            end
%         end
%         
       if length(ind_lines_intersect) == 2
            valency = 2;
            return;
       end
           
           
       valency = estimateValenceLinesGroup(ind_lines_intersect, collinear_strokes_pairs,  intersections.strokes_indices, strokes_topology);
        
%        if valency >=3
%                 for i = ind_lines_intersect
%                     plot(cat(1,strokes_topology(i).points2D(:).x), cat(1,strokes_topology(i).points2D(:).y));
%                 end
%            end
               
%         plot(cat(1,strokes_topology(i).points2D(:).x), cat(1,strokes_topology(i).points2D(:).y));
%         for j = ind_lines_intersect
%             plot(cat(1,strokes_topology(j).points2D(:).x), cat(1,strokes_topology(j).points2D(:).y), ':', 'LineWidth', 2);
%         end
%         plot(intersections.coordinates2D(int_ind, 1),...
%                 intersections.coordinates2D(int_ind, 2), '*b');  
%
end




function valency = estimateValenceLinesGroup(ind_lines_intersect, collinear_strokes_pairs, intersections_pairs, strokes_topology)
       
       combos = nchoosek(1:length(ind_lines_intersect),2);
       combos = orderPairs(combos);
       pairs = ind_lines_intersect(combos);
%        pairs = orderPairs(pairs);
       intersections_pairs = orderPairs(intersections_pairs);
       

       ind_compute_collinearity = find(~ismember(pairs, intersections_pairs, 'rows'));

       for pairs_ind = 1:length(ind_compute_collinearity)
            ind_s1 = pairs(ind_compute_collinearity(pairs_ind),1);
            ind_s2 = pairs(ind_compute_collinearity(pairs_ind),2);
    
            if evaluateIfStrokesCollinear(strokes_topology, ind_s1, ind_s2)
%                 disp([ind_s1, ind_s2]);
                collinear_strokes_pairs(end+1,:) = sort([ind_s1, ind_s2]);
            end

       end
       
       %        [sb_xy, se_xy, dir] = convertSegIndToSegInterval(stroke, seg_num);
       
       
       mask_edges = ismember(pairs, collinear_strokes_pairs, 'rows');

       
%        mask_edges = (ismember(ind_lines_intersect, collinear_strokes_pairs, 'rows'));
       
%        old_vals = unique(ind_lines_intersect);
       new_vals = 1:length(ind_lines_intersect);
% %        
%        [~, pairs] = ismember(ind_lines_intersect,old_vals);
%        
%        
%        s = pairs(mask_edges,1);
%        t = pairs(mask_edges,2);
       
       
       s = combos(mask_edges,1);
       t = combos(mask_edges,2);
       
%        try
       G = graph([s' new_vals],[t' new_vals ]); 
       
        
%        catch
%            disp('');
%        end
       valency = length(unique(conncomp(G)));
end