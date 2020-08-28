lists;

%% Read json data for original strokes:
% keep only starigth stroke and those that have imformation about the
% assignement.
path_json = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\full_json';
percentage_low_score_dave_14 = [];
percentage_low_score_jerry_14 = [];

for i = 1:size(list_designers_paper_14,1)
    desinger{i} = list_designers_paper_14{i,1};
    object{i} = list_designers_paper_14{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough_full.json'];
    
   [strokes_topology_rough, ~, ~] = ...
            readReconstructionJson( fullfile(path_json, filename)); 
   
    % dave
    
   set_removed = strokes_numbers_dave_removed{i};
   percentage_low_score_dave_14(i) = getConfidenceScoreVals(set_removed, strokes_topology_rough);
   
   
   set_removed = strokes_numbers_jerry_removed{i};
   percentage_low_score_jerry_14(i) = getConfidenceScoreVals(set_removed, strokes_topology_rough);
   
end

fprintf('Mean_percentage_low_score_dave = %.2f\n', mean(percentage_low_score_dave_14 (~isnan(percentage_low_score_dave_14 ))));
fprintf('Median_percentage_low_score_dave = %.2f\n', median(percentage_low_score_dave_14 (~isnan(percentage_low_score_dave_14 ))));

fprintf('Mean_percentage_low_score_jerry = %.2f\n', mean(percentage_low_score_jerry_14(~isnan(percentage_low_score_jerry_14))));
fprintf('Median_percentage_low_score_jerry = %.2f\n', median(percentage_low_score_jerry_14(~isnan(percentage_low_score_jerry_14))));

maks_average = find(~isnan(percentage_low_score_dave_14 ) & ~isnan(percentage_low_score_jerry_14));
maks_dave = find(~isnan(percentage_low_score_dave_14 ) & isnan(percentage_low_score_jerry_14));
maks_jerry = find(isnan(percentage_low_score_dave_14 ) & ~isnan(percentage_low_score_jerry_14));

percentage_low_score_14 = NaN*ones(size(percentage_low_score_dave_14 ));
percentage_low_score_14(maks_average) = 0.5*(percentage_low_score_dave_14 (maks_average)+percentage_low_score_jerry_14(maks_average));
percentage_low_score_14(maks_dave) = percentage_low_score_dave_14 (maks_dave);
percentage_low_score_14(maks_jerry) = percentage_low_score_jerry_14(maks_jerry);

fprintf('Mean_percentage_low_score_14 = %.2f\n', mean(percentage_low_score_14(~isnan(percentage_low_score_14))));
fprintf('Median_percentage_low_score_14 = %.2f\n', median(percentage_low_score_14(~isnan(percentage_low_score_14))));


%% supp
percentage_low_score_yuan_supp = [];
percentage_low_score_jerry_supp = [];

for i = 1:size(list_designers_supplemental,1)
    desinger{i} = list_designers_supplemental{i,1};
    object{i} = list_designers_supplemental{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough_full.json'];
    
   [strokes_topology_rough, ~, ~] = ...
            readReconstructionJson( fullfile(path_json, filename)); 
   
    % yuan
    
   set_removed = strokes_numbers_yuan_removed_supp{i};
    
   pcl = getConfidenceScoreVals(set_removed, strokes_topology_rough);
   percentage_low_score_yuan_supp(i) = pcl;
   
   % jerry
   set_removed = strokes_numbers_jerry_removed_supp{i};
   
   pcl = getConfidenceScoreVals(set_removed, strokes_topology_rough);
   percentage_low_score_jerry_supp(i) = pcl;
end

fprintf('Mean_percentage_low_score_yuan_supp = %.2f\n', mean(percentage_low_score_yuan_supp(~isnan(percentage_low_score_yuan_supp))));
fprintf('Median_percentage_low_score_yuan_supp = %.2f\n', median(percentage_low_score_yuan_supp(~isnan(percentage_low_score_yuan_supp))));

fprintf('Mean_percentage_low_score_jerry_supp = %.2f\n', mean(percentage_low_score_jerry_supp(~isnan(percentage_low_score_jerry_supp))));
fprintf('Median_percentage_low_score_jerry_supp = %.2f\n', median(percentage_low_score_jerry_supp(~isnan(percentage_low_score_jerry_supp))));

maks_average = find(~isnan(percentage_low_score_yuan_supp) & ~isnan(percentage_low_score_jerry_supp));
maks_dave = find(~isnan(percentage_low_score_yuan_supp) & isnan(percentage_low_score_jerry_supp));
maks_jerry = find(isnan(percentage_low_score_yuan_supp) & ~isnan(percentage_low_score_jerry_supp));

percentage_low_score_supp = NaN*ones(size(percentage_low_score_yuan_supp));
percentage_low_score_supp(maks_average) = 0.5*(percentage_low_score_yuan_supp(maks_average)+percentage_low_score_jerry_supp(maks_average));
percentage_low_score_supp(maks_dave) = percentage_low_score_yuan_supp(maks_dave);
percentage_low_score_supp(maks_jerry) = percentage_low_score_jerry_supp(maks_jerry);


v_supp = percentage_low_score_supp(~isnan(percentage_low_score_supp));
v14 = percentage_low_score_14(~isnan(percentage_low_score_14));

fprintf('Mean_percentage_low_score_supp = %.2f\n', mean(v_supp));
fprintf('Median_percentage_low_score_supp = %.2f\n', median(percentage_low_score_supp(~isnan(percentage_low_score_supp))));

fprintf('Mean_percentage_low_score_all = %.2f\n', mean( [v_supp  v14]));
fprintf('Median_percentage_low_score_all = %.2f\n', median([v_supp  v14]  ));
