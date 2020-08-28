function [agreement, solution_labeling, solution_our] = compareSolutions(solution_labeling, agreement, inds_new, intersections_rec,...
                        folder_save,sketch_height,strokes)

 
 mask_pair = ~isnan(inds_new);
 solution_labeling = solution_labeling(mask_pair);
 agreement = agreement(mask_pair);
 inds_marked = find(mask_pair);
 inds_new = inds_new(mask_pair);

 mask_likely = cat(1,cat(1,intersections_rec(inds_new).likely));
 inds_marked = inds_marked(mask_likely);
 inds_new = inds_new(mask_likely);
 solution_labeling = solution_labeling(mask_likely);
 agreement = agreement(mask_likely);
 
 
 solution_our = cat(1,intersections_rec(inds_new).is_active);
 
 
 mask_nonempty = arrayfun(@(x) ~isempty(x.is_active),intersections_rec(inds_new));
 
 inds_new = inds_new(mask_nonempty);
 
 solution_labeling = solution_labeling(mask_nonempty);
 agreement = agreement(mask_nonempty);
 
 intersections_coordinates = cat(1,intersections_rec(inds_new).coordinates2D);
 
 plot(intersections_coordinates(:,1), intersections_coordinates(:,2), 'r*')
 
 mask_agreement_75 =  agreement > 0.75;
 mask_agreement_50 =  (agreement > 0.5) & agreement < 0.75;
 
 percentage_same = sum((solution_our == solution_labeling'))./length(solution_labeling);
 
 percentage_same_75 = sum((solution_our(mask_agreement_75) == ...
                        solution_labeling(mask_agreement_75)'))./...
                        length(solution_labeling(mask_agreement_75));
                    
 percentage_same_50 = sum((solution_our(mask_agreement_50) == ...
                        solution_labeling(mask_agreement_50)'))./...
                        length(solution_labeling(mask_agreement_50));    
                    
                    
 fprintf('full = %.3f, 75 = %.3f, 50 = %.3f\n', ...
            percentage_same, ...
            percentage_same_75, ...
            percentage_same_50)                    
 
% save_as_svg_our_vs_human(strokes,intersections_coordinates,solution_labeling, solution_our, folder_save, sketch_height, 'human_vs_our.svg')
% save_as_svg_ours(strokes,intersections_coordinates,solution_our, folder_save, sketch_height, 'our.svg')

end